#include "../src/Grid/Grid"
#include "../src/Utility/hashmap.h"

using namespace Larus::Grid;

inline void test_show_node() {
	typedef Node<CooValueType, ValueType, 3> Node;

	Node* pnode = new Node(NULL_PTR, 0, 0, 0, 0, //
			0.5, 0.5, //
			0.5, 0.5, //
			0.5, 0.5);
	//

	vtk_show(pnode);

	delete pnode;

}

inline void test_show_grid() {
	typedef Grid<CooValueType, ValueType, 3> Grid;
	Grid g( 3, 0, 1, //
			3, 0, 1, //
			3, 0, 1);
	//
	vtk_show(g);
}


inline void test_hashmap(){
	Larus::HashMap<double, double> hm;
	for(int i = 0; i<10; ++i){
		hm.Insert(i, 500+i);
	}
	for(auto iter=hm.begin(); iter!=hm.end(); ++iter){
		std::cout<< (*iter) <<"\n";
	}
	std::cout<< hm<<"\n";
	std::cout<< hm.Get(1) <<"\n";
}
