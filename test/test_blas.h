#include "../src/IO/IO_gnuplot.h"
#include "../src/IO/IO_vtk.h"
#include "../src/Geometry/Triangle.h"
#include "../src/Geometry/Line.h"
#include "../src/Geometry/Polygon.h"
#include "../src/Utility/Array.h"
#include "../src/Geometry/Relation.h"
#include "../src/Geometry/Plane.h"
#include "../src/Grid/SPTreeNode.h"
#include "../src/Grid/SPTree.h"
#include "../src/Grid/Cell.h"
#include "../src/Grid/Forest.h"
#include "../src/Algebra/Space.h"
#include "../src/Calculation/VOF.h"
#include "../src/Calculation/Levelset.h"
#include "../src/Calculation/Scalar.h"
#include "../src/Algebra/Arithmetic.h"
#include "../src/Algebra/BLAS_Level1.h"
#include "../src/Algebra/Cube.h"
#include "../src/TypeDef.h"
#include "../src/Algebra/Interpolation.h"
#include "../src/Algebra/Expression.h"
#include "../src/Grid/Dimension.h"
#include "../src/IO/IO_gerris.h"
#include <list>
#include <math.h>
#include <sstream>

using namespace std;
namespace Larus
{
void test_blas_1()
{
	cout << "test blas ==================\n";
	int n = 7;
	arrayList ax(n);
	ax.assign(3);
	arrayList ay(n);
	ay.assign(5);
	ax.show();
	ay.show();
	cout << "==================\n";
	cout << nrmp(ax, 0.5) << " dot \n";
	ax.show();
	ay.show();
}

void test_simd()
{

	int n = 100000000;
	arrayListV<float> af(n);
	af.assign(1);
	arrayListV<float> af2(n);
	af2.assign(2);
	tick_t tb = Clock::Tick();
	af2 = af2 + af;
	//for (int i = 0; i < af2.size(); i++)
	//{
	//	cout << af2[i] << endl;
	//}
	tick_t te = Clock::Tick()-tb;
	cout << "Time float   " << Clock::TicksToMillisecondsF(te)
			<< endl;
	cout << "=================================" << endl;
	//Clock::Sleep(100);

	arrayListV<double> aff(n);
	aff.assign(4);
	arrayListV<double> aff2(n);
	aff2.assign(5);
	tick_t tbf = Clock::Tick();
	aff2 = aff2 + aff;
	//for (int i = 0; i < aff2.size(); i++)
	//{
	//	cout << aff2[i] << endl;
	//}
	tick_t tef = Clock::Tick() - tbf;
	cout << "Time double  " << Clock::TicksToMillisecondsF(tef)
			<< endl;
}
}
