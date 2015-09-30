#ifndef DATA_H_
#define DATA_H_

#include "../TypeDef.h"
#include "grid_def.h"
#include "../Utility/ArrayList.h"

namespace Larus {
namespace Grid {
template<typename VALUE, size_t DIM>
class Data {
public:
	static const size_t Dim = DIM;
	static const size_t NumFaces = DIM + DIM;
	static const size_t NumVertexes = (DIM == 3) ? 8 : (DIM + DIM); //

	typedef VALUE value_t;
	typedef Data<VALUE, DIM> self;

	typedef void (*pfunction)(self*, utPointer);

protected:
	arrayListV<value_t> _center;
	arrayListV<value_t> _face[NumFaces];
	arrayListV<value_t> _vertex[NumVertexes];
	utPointer untype;
public:
	Data() {
		untype = NULL_PTR;
	}
	Data(const size_t& nc, const size_t& nf, const size_t& nv, utPointer utp) :
			_center(nc) {
		for (int i = 0; i < NumFaces; ++i) {
			_face[i].reconstruct(nf);
		}
		for (int i = 0; i < NumVertexes; ++i) {
			_face[i].reconstruct(nv);
		}
		untype = utp;
	}

	bool is_empty() const {
		bool res = true;
		res = res && (_center.size() == 0);
		for (int i = 0; i < NumFaces; ++i) {
			res = res && (_face[i].size() == 0);
		}
		for (int i = 0; i < NumVertexes; ++i) {
			res = res && (_vertex[i].size() == 0);
		}
		return res;
	}

	void show_info() const {
		std::cout << "center data:" << this->_center.size() << "\n";
		std::cout << "face data   :" << "\n";
		for (int i = 0; i < NumFaces; ++i) {
			std::cout << "    " << i << "       :" << _face[i].size()
					<< "\n";
		}
		for (int i = 0; i < NumVertexes; ++i) {
			std::cout << "    " << i << "       :" << _vertex[i].size()
					<< "\n";
		}
	}

};
}
}

#endif /* DATA_H_ */
