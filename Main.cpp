#include <Siv3D.hpp> // OpenSiv3D v0.4.0

using Index = size_t;
using Vertices = Array<Vec2>;

class QuadTreeVertices
{
public:
	enum Dir { TL = 0, TR, BR, BL };

	QuadTreeVertices() = default;

	QuadTreeVertices(const RectF& scope)
		: scope(scope)
		, depth(0)
	{}

	void init(const RectF& currentScope, int currentDepth = 0)
	{
		depth = currentDepth;
		scope = currentScope;
		children.clear();
		is.clear();
		vs.clear();
	}

	void add(const Vertices& vertices, Index index)
	{
		const auto addToChild = [&](Index i)
		{
			const Vec2& vertex =  vertices[i];
			const Vec2 center = scope.center();
			//Left
			if (vertex.x < center.x)
			{
				children[vertex.y < center.y ? TL : BL].add(vertices, index);
			}
			//Right
			else
			{
				children[vertex.y < center.y ? TR : BR].add(vertices, index);
			}
		};

		//既に分割済みの場合
		if (!children.empty())
		{
			addToChild(index);
		}
		//既に頂点を1つ以上持っていて最大深度より低ければ子要素に全て渡す
		else if (!is.empty() && depth < MaxDepth)
		{
			const Vec2 childWidth = scope.size * 0.5;
			const int childDepth = depth + 1;

			children.resize(4);
			children[TL].init(RectF(scope.pos + Vec2(0, 0) * childWidth, childWidth), childDepth);
			children[TR].init(RectF(scope.pos + Vec2(1, 0) * childWidth, childWidth), childDepth);
			children[BR].init(RectF(scope.pos + Vec2(1, 1) * childWidth, childWidth), childDepth);
			children[BL].init(RectF(scope.pos + Vec2(0, 1) * childWidth, childWidth), childDepth);

			for (auto i : is)
			{
				addToChild(i);
			}
			is.clear();
			vs.clear();

			addToChild(index);
		}
		//空のノードであれば頂点を1つ持つリーフノードになる
		else
		{
			is.push_back(index);
			vs.push_back(vertices[index]);
		}
	}

	static int MaxDepth;

	Optional<Vec2> centroid()const
	{
		if (countOfChildren == 0)
		{
			return none;
		}
		return sumOfChildren / countOfChildren;
	}

	const Vec2& size()const
	{
		return scope.size;
	}

	bool isTree()const
	{
		return !isLeaf();
	}

	bool isLeaf()const
	{
		return children.empty();
	}

	void update()
	{
		updateImpl();
	}

	Vec2 repulsiveForce(const Vec2& xi, Index index, const std::function<Vec2(const Vec2&, const Vec2&)>& fr)const
	{
		const size_t S = countOfChildren;
		if (1 <= S)
		{
			const Vec2 xS = sumOfChildren / S;

			//Barnes-Hut opening criterion
			const double theta = 1.2;
			const double dS = Max(scope.w, scope.h);
			const bool usableAsSuperNode = ((dS * dS) / (xi - xS).lengthSq() <= theta * theta);

			if (usableAsSuperNode)
			{
				return S * fr(xi, xS);
			}
			else if (isTree())
			{
				Vec2 f = Vec2::Zero();
				for (auto& child : children)
				{
					f += child.repulsiveForce(xi, index, fr);
				}
				return f;
			}
			else
			{
				Vec2 f = Vec2::Zero();
				for (auto [i, v] : Indexed(vs))
				{
					if (is[i] == index)
					{
						continue;
					}

					f += fr(xi, v);
				}
				return f;
			}
		}

		return Vec2::Zero();
	}

private:
	std::pair<Vec2, size_t> updateImpl()
	{
		sumOfChildren = Vec2::Zero();
		countOfChildren = 0;

		for (auto& child : children)
		{
			auto [v, s] = child.updateImpl();
			sumOfChildren += v;
			countOfChildren += s;
		}

		if (!vs.empty())
		{
			sumOfChildren += vs.reduce1([](const Vec2& a, const Vec2& b) { return a + b; });
			countOfChildren += vs.size();
		}

		return { sumOfChildren, countOfChildren };
	}

	RectF scope;
	Array<QuadTreeVertices> children;
	Array<Vec2> vs;
	Array<Index> is;
	int depth;

	Vec2 sumOfChildren;
	size_t countOfChildren;
};

int QuadTreeVertices::MaxDepth;

class Graph
{
public:
	using Indices = std::unordered_set<Index>;
	using Links = std::vector<Indices>;

	Graph() = default;

	//.mtxフォーマットを開く
	Graph(const FilePath& path)
	{
		TextReader reader(path);
		if (!reader.isOpened())
		{
			Logger << U"エラー：ファイル \"" << path << U"\" が開けません。";
			return;
		}

		size_t nodeSize = 0;
		size_t linkSize = 0;
		while (auto lineOpt = reader.readLine())
		{
			if (const bool isComment = lineOpt.value().starts_with(U"%%"))
			{
				continue;
			}

			const auto nums = lineOpt.value().split(U' ').map([](const String& str) { return ParseInt<int>(str); });
			if (nums.size() == 3)
			{
				if (const bool notSquareMatrix = (nums[0] != nums[1]))
				{
					Logger << U"エラー：グラフのフォーマットが不正です。";
					return;
				}

				nodeSize = nums[0];
				linkSize = nums[2];
				break;
			}
		}

		links.resize(nodeSize);
		while (auto lineOpt = reader.readLine())
		{
			if (const bool isComment = lineOpt.value().starts_with(U"%%"))
			{
				continue;
			}

			auto nums = lineOpt.value().split(U' ').map([](const String& str) { return ParseInt<int>(str); });
			if (nums.size() == 2)
			{
				const Index indexA = nums[0] - 1;
				const Index indexB = nums[1] - 1;
				links[indexA].insert(indexB);
				links[indexB].insert(indexA);
			}
		}

		nodes.resize(nodeSize);
		reset();
	}

	void reset()
	{
		const double initialRadius = 2500.0;
		for (auto& node : nodes)
		{
			node = RandomVec2(Circle(Vec2::Zero(), initialRadius));
		}

		timeStep = 1.0;
		energy = DBL_MAX;
		progress = 0;
		maxTreeDepth = 8;

		//model parameter
		C = 0.2;
		K = 100.0;
	}

	void update()
	{
		Stopwatch watch(true);
		for (int i = 0; i < 7; ++i)
		{
			updateWithBarnesHutForce();
		}
		//Console << watch.msF();
	}

	void draw()const
	{
		for (auto [i, linksTo] : Indexed(links))
		{
			const Vec2 p0 = nodes[i];
			for (auto linkTo : linksTo)
			{
				const Vec2 p1 = nodes[linkTo];
				Line(p0, p1).draw(2,Palette::White);
			}
		}
	}

	Vec2 centroid()const
	{
		return nodes.reduce([](const Vec2& value, const Vec2& pos) { return value + pos; }, Vec2::Zero()) / nodes.size();
	}

	RectF aabb()const
	{
		Vec2 minPos(DBL_MAX, DBL_MAX), maxPos(-DBL_MAX, -DBL_MAX);
		nodes.each([&](const Vec2& p)
			{
				minPos.x = Min(minPos.x, p.x);
				minPos.y = Min(minPos.y, p.y);
				maxPos.x = Max(maxPos.x, p.x);
				maxPos.y = Max(maxPos.y, p.y);
			});

		return RectF(minPos, maxPos - minPos);
	}

private:
	void updateWithBarnesHutForce()
	{
		double energy0 = energy;
		energy = 0;

		const double eps = 1.e-4;
		quadTree.init(aabb().stretched(eps));

		QuadTreeVertices::MaxDepth = maxTreeDepth;
		for (auto i : step(nodes.size()))
		{
			quadTree.add(nodes, i);
		}

		quadTree.update();
		
		for (auto [i, linksTo] : Indexed(links))
		{
			Vec2& p0 = nodes[i];

			Vec2 force = Vec2::Zero();

			const double p = 2.0;
			const double coeff = -C * pow(K, 1.0 + p);
			const auto fr = [=](const Vec2& xi, const Vec2& xj)
			{
				const Vec2 v = (xj - xi);
				return coeff * v / pow(v.lengthSq(), (1.0 + p) * 0.5);
			};

			// 反発力の計算: O(nlogn)
			force += quadTree.repulsiveForce(p0, i, fr);

			/*
			// 反発力の計算: O(|V|^2)
			for (Index j = 0; j < nodes.size(); ++j)
			{
				if (i == j)
				{
					continue;
				}

				force += fr(p0, nodes[j]);
			}
			*/

			//引力の計算
			for (auto linkTo : linksTo)
			{
				if (i == linkTo)
				{
					continue;
				}

				const Vec2 p1 = nodes[linkTo];
				const Vec2 v = (p1 - p0);
				force += v.length() * v / K;
			}

			p0 += timeStep * force.normalized();
			energy += force.lengthSq();
		}
	}

	Links links;
	Vertices nodes;
	double timeStep;
	double energy;
	int progress;
	int maxTreeDepth;
	double C;
	double K;

	QuadTreeVertices quadTree;
};

void Main()
{
	int windowWidth = 1000;
	Window::Resize(windowWidth, windowWidth);

	//download from: https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/jagmesh/jagmesh1.html
	Graph graph(U"jagmesh1.mtx");

	Camera2D camera;
	camera.setTargetScale(windowWidth / 5200.0);
	
	while (System::Update())
	{
		if (KeySpace.down())
		{
			graph.reset();
		}

		{
			auto t = camera.createTransformer();

			graph.update();
			graph.draw();
		}

		const Vec2 centroid = graph.centroid();
		/*
		const RectF aabb = graph.aabb();

		const double desiredWidth = 1.1 * 2.0 * Max({
			abs(aabb.tl().x - centroid.x),
			abs(aabb.tl().y - centroid.y),
			abs(aabb.br().x - centroid.x),
			abs(aabb.br().y - centroid.y)
			});

		camera.setTargetScale(windowWidth / desiredWidth);
		*/
		camera.setTargetCenter(centroid);

		camera.update();
		camera.draw();
	}
}