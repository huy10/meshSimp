// MeshSimp.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "SimpleObject.h"
#include "Vec3f.h"

#include "MeshSimp.h"
#include <algorithm>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/vector.hpp> 
#include <boost/numeric/ublas/io.hpp> 
#include <boost/numeric/ublas/vector_proxy.hpp> 
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/triangular.hpp> 
#include <boost/numeric/ublas/lu.hpp> 
 

using namespace std;
using namespace SimpleOBJ;
// compute Q
void computeQ(int index,  int m_nVertices, int m_nTriangles, Vec3f * v, Array<int,3>* Triangles, boost::numeric::ublas::matrix <float> &Q)
{
	float a = 0, b = 0, c = 0, d = 0;
	for (unsigned i = 0; i < Q.size1 (); ++ i)
        for (unsigned j = 0; j < Q.size2 (); ++ j)
            Q(i, j) = 0;

	for (int i = 0; i < m_nTriangles ; i++)
	{
		if (Triangles[i][0] == index || Triangles[i][1] == index || Triangles[i][2] == index)
		{
			
			boost::numeric::ublas::matrix <float> K(4,4);      // this is Kp
			float x_0 = v[Triangles[i][0]].x;
			float y_0 = v[Triangles[i][0]].y;
			float z_0 = v[Triangles[i][0]].z;
			float x_1 = v[Triangles[i][1]].x;
			float y_1 = v[Triangles[i][1]].y;
			float z_1 = v[Triangles[i][1]].z;
			float x_2 = v[Triangles[i][2]].x;
			float y_2 = v[Triangles[i][2]].y;
			float z_2 = v[Triangles[i][2]].z;
			if (   z_2*(y_0-y_1) != y_2*(z_0-z_1) )
			{
					float A = - ((x_0-x_1)- x_2*(z_0-z_1)/z_2) /  ((y_0-y_1) - y_2*(z_0-z_1)/z_2);
					float B = - (x_2 + y_2 * A) / z_2;
					a = sqrt(1/( 1 + pow(A,2) + pow(B,2)));
					b = A * a;
					c = B * a;
					d = - (a*x_0 + b*y_0 + c*z_0);
			}
			else
			{
				a = 0; 
				b = 1;     //z_0 / sqrt(pow(z_0,2) + pow(y_0,2));
				c = 0 ;     // / sqrt(pow(z_0,2) + pow(y_0,2));
				d = 0;
			}
			
			float p [4] = {a,b,c,d};
			for (unsigned i = 0; i < Q.size1 (); ++ i)
				 for (unsigned j = 0; j < Q.size2 (); ++ j)
					K(i, j) = p [i] * p [j];
			//cout << K << endl;
			Q += K;
		}
	}

	//std::cout << "index = " << index << endl;
	//cout << Q << endl;
}

// this is the inverse function
namespace ublas = boost::numeric::ublas;
 
// template<class E1, class E2>
//InvertMatrix inverse (matrix_expression<E1> &e1, matrix_expression<E2> &e2) {
//    typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
//    typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
//    typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;
//    BOOST_UBLAS_CHECK (e1 ().size1 () == e2 ().size1 (), bad_size ());
//    BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size2 (), bad_size ());
//    size_type size = e1 ().size1 ();
//    for (size_type n = 0; n < size; ++ n) {
//       // processing column n
//       // find the row that has the largest number at this column (in absolute value)
//       size_type best_row = index_norm_inf(row(e1(), n));
//       // check wether this number is'nt zero
//         BOOST_UBLAS_CHECK (e1 () (best_row, n) != value_type (), singular ());
//       { // swap this row with the n-th row
//          vector<value_type> temp = row(e1(), best_row);
//          row(e1(), n) = row(e1(), best_row);
//          row(e1(), best_row) = temp;
//       }
//       // do the same on the destination matrix
//       { // swap this row with the n-th row
//          vector<value_type> temp = row(e2(), best_row);
//          row(e2(), n) = row(e2(), best_row);
//          row(e2(), best_row) = temp;
//       }
//       // now eliminate all elements below and above this row
//       for (size_type i = 0; i < size; ++ i)
//          if (i!=n) {
//             value_type t = -e1 () (i, n)/ e1 () (n, n);
//             row(e1(), i) += t*row(e1(), n);
//             row(e2(), i) += t*row(e2(), n);
//          } else {
//             value_type t = 1 / e1 () (i, n);
//             row(e1(), i) *= t;
//             row(e2(), i) *= t;
//          }
//    }
// }

 //template<class E1, class E2>
 //void inverse (matrix_expression<E1> &e1, matrix_expression<E2> &e2) {
 //   typedef BOOST_UBLAS_TYPENAME E2::size_type size_type;
 //   typedef BOOST_UBLAS_TYPENAME E2::difference_type difference_type;
 //   typedef BOOST_UBLAS_TYPENAME E2::value_type value_type;
 //   BOOST_UBLAS_CHECK (e1 ().size1 () == e2 ().size1 (), bad_size ());
 //   BOOST_UBLAS_CHECK (e1 ().size2 () == e2 ().size2 (), bad_size ());
 //   size_type size = e1 ().size1 ();
 //   for (size_type n = 0; n < size; ++ n) {
 //      // processing column n
 //      // find the row that has the largest number at this column (in absolute value)
 //      size_type best_row = index_norm_inf(row(e1(), n));
 //      // check wether this number is'nt zero
 //        BOOST_UBLAS_CHECK (e1 () (best_row, n) != value_type (), singular ());
 //      { // swap this row with the n-th row
 //         vector<value_type> temp = row(e1(), best_row);
 //         row(e1(), n) = row(e1(), best_row);
 //         row(e1(), best_row) = temp;
 //      }
 //      // do the same on the destination matrix
 //      { // swap this row with the n-th row
 //         vector<value_type> temp = row(e2(), best_row);
 //         row(e2(), n) = row(e2(), best_row);
 //         row(e2(), best_row) = temp;
 //      }
 //      // now eliminate all elements below and above this row
 //      for (size_type i = 0; i < size; ++ i)
 //         if (i!=n) {
 //            value_type t = -e1 () (i, n)/ e1 () (n, n);
 //            row(e1(), i) += t*row(e1(), n);
 //            row(e2(), i) += t*row(e2(), n);
 //         } else {
 //            value_type t = 1 / e1 () (i, n);
 //            row(e1(), i) *= t;
 //            row(e2(), i) *= t;
 //         }
 //   }
 //}
 /* Matrix inversion routine.
    Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
 //template<class T>
 bool InvertMatrix (const ublas::matrix<float >& input, ublas::matrix<float>& inverse) 
 {
 	using namespace boost::numeric::ublas;
 	typedef permutation_matrix<std::size_t> pmatrix;
 	// create a working copy of the input
 	matrix<float> A(input);
 	// create a permutation matrix for the LU-factorization
 	pmatrix pm(A.size1());
 
 	// perform LU-factorization
 	int res = lu_factorize( A, pm );
        if( res != 0 ) return false;
 
 	// create identity matrix of "inverse"
 	inverse.assign(ublas::identity_matrix<float>(A.size1()));

 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);
 
 	return true;
 }

 bool compareCost(_pair a, _pair b)
 {
	 return a.cost < b.cost ; 
 }

 bool isEdge(int indexa, int indexb, Array<int,3>* Triangles, int m_nTriangles)
 {
	 for (unsigned int i = 0; i < m_nTriangles; i++)
	 {
		 if (Triangles[i][0] == indexa || Triangles[i][1] == indexa || Triangles[i][2] == indexa)
			 if (Triangles[i][0] == indexb || Triangles[i][1] == indexb|| Triangles[i][2] == indexb)
				{
					//cout << Triangles[i][0]+1 << "    "<< Triangles[i][1] +1<<  "  "<< Triangles[i][2]+1<<endl;  
					return true;
			 }
		 
	 }
	 return false;
 }
 bool m_equi( _pair a, _pair b)
{
	if (a.v1 ==  b.v1 && a.v2 == b.v2)
		return true;
	else if (a.v1 ==  b.v2 && a.v2 == b.v1)
		return true;
	else 
		return false;
}
	
int _tmain(int argc, _TCHAR* argv[])
{
		
	

	SimpleOBJ::CSimpleObject mOBJ;

	bool loadSuccess =  mOBJ.LoadFromObj(fileName);


	//mOBJ.m_pVertexList[0].
	// compute Q for all vertexes
	
	boost::numeric::ublas::matrix <float> * Qs = new boost::numeric::ublas::matrix <float> [mOBJ.m_nVertices]; 


	for (int i = 0; i < mOBJ.getn_v(); i++)
	{
		float a = 0, b = 0, c = 0, d = 0;
		boost::numeric::ublas::matrix <float> Q (4,4);
		computeQ( i,    mOBJ.m_nVertices,      mOBJ.m_nTriangles,    mOBJ.m_pVertexList,     mOBJ.m_pTriangleList,   Q);
		Qs[i] = Q;
		//cout << Qs[i] << endl;
		
	}


	// select valid pairs
	vector<_pair> pairs;      // the index of v1, v2, float cost and Vec3f vbar

	for (unsigned int i = 0; i < mOBJ.m_nTriangles; i++)
	{
		_pair edge;
		edge.v1 = mOBJ.m_pTriangleList[i][0];
		edge.v2 = mOBJ.m_pTriangleList[i][1];
	
		if (edge.v1 > edge.v2)
		{
			unsigned int temp = edge.v1;
			edge.v1 = edge.v2;
			edge.v2 = temp;
		}
		/*bool exist = false;
		for (unsigned j = 0; j < pairs.size(); j++)
		{
			if (edge.v1 == pairs.at(j).v1 && edge.v2 == pairs.at(j).v2 )
				exist = true;
		}
		if (exist == false)*/
			pairs.push_back(edge);

		edge.v1 = mOBJ.m_pTriangleList[i][0];
		edge.v2 = mOBJ.m_pTriangleList[i][2];
		if (edge.v1 > edge.v2)
		{
			unsigned int temp = edge.v1;
			edge.v1 = edge.v2;
			edge.v2 = temp;
		}
		/*exist = false;
		for (unsigned j = 0; j < pairs.size(); j++)
		{
			if (edge.v1 == pairs.at(j).v1 && edge.v2 == pairs.at(j).v2 )
				exist = true;
		}
		if (exist == false)*/
			pairs.push_back(edge);

		edge.v1 = mOBJ.m_pTriangleList[i][1];
		edge.v2 = mOBJ.m_pTriangleList[i][2];
		if (edge.v1 > edge.v2)
		{
			unsigned int temp = edge.v1;
			edge.v1 = edge.v2;
			edge.v2 = temp;
		}
		/*exist = false;
		for (unsigned j = 0; j < pairs.size(); j++)
		{
			if (edge.v1 == pairs.at(j).v1 && edge.v2 == pairs.at(j).v2 )
				exist = true;
		}
		if (exist == false)*/
			pairs.push_back(edge);
	}


#ifdef Debug
	cout << "pairs.size() = " << pairs.size() << endl;
	cout << "edge push back finished" << endl;
#endif

	for (unsigned int i = 0; i < mOBJ.m_nVertices; i++)
		for (unsigned int j = i+1; j < mOBJ.m_nVertices; j++)
						
		{
			
			Vec3f minusVec = mOBJ.m_pVertexList[i] - mOBJ.m_pVertexList[j];
			float L2dis  = minusVec.L2Norm_Sqr();
			float t = 0.0;
			if ( L2dis < t)
				{
					_pair edge;
					edge.v1 = i;
					edge.v2 = j;
					if (edge.v1 > edge.v2)
					{
						unsigned int temp = edge.v1;
						edge.v1 = edge.v2;
						edge.v2 = temp;
					}
					pairs.push_back(edge);
				}
			
		}
		
	std::vector<_pair>::iterator it;
	cout << "pairs.size() = " << pairs.size() << endl;
	it = unique(pairs.begin(),pairs.end() ,  m_equi);
	pairs.resize( std::distance(pairs.begin(),it) );
#ifdef Debug
	//for (unsigned int i = 0; i < pairs.size();i++)
	//{
	//	cout << "v1:" << pairs.at(i).v1 << "  v2:" << pairs.at(i).v2 << endl;
	//}

	cout << "pairs.size() = " << pairs.size() << endl;
	cout << "Pair selection finished" << endl;
#endif
	// now we have got all valid pairs
	for (unsigned int i=0; i < pairs.size();i++)
	{
		// compute cost  = Q*vbar
		_pair edge = pairs.at(i);
		boost::numeric::ublas::matrix <float> Q1 (4,4);
		Q1 = Qs[edge.v1];
		
		assert(edge.v1 < edge.v2);
		boost::numeric::ublas::matrix <float> Q2 (4,4);
		Q2 = Qs[edge.v2];
		boost::numeric::ublas::matrix <float> Q (4,4);
		Q = Q1 + Q2;
				
		boost::numeric::ublas::matrix <float> Q_inv (4,4);
		bool invertable = InvertMatrix(Q,Q_inv);
		
		if (invertable == true)
		{
			
			pairs.at(i).vbar.push_back(Q_inv(0,3) / Q_inv(3,3));
			pairs.at(i).vbar.push_back(Q_inv(1,3) / Q_inv(3,3));
			pairs.at(i).vbar.push_back(Q_inv(2,3) / Q_inv(3,3));
			pairs.at(i).vbar.push_back(1);
			

			boost::numeric::ublas::vector <float> vBar (4);
			vBar(0) = pairs.at(i).vbar.at(0);
			vBar(0) = pairs.at(i).vbar.at(1);
			vBar(0) = pairs.at(i).vbar.at(2);
			vBar(0) = pairs.at(i).vbar.at(3);

			boost::numeric::ublas::vector <float> vBart (4);
			vBart  = prod(Q, vBar);
			//cout << vBar(0) << "  " << vBar(1) << "   "<<vBar(2) << "  " << vBar(3) << endl;
			float cost = 0;
			for (unsigned int i = 0; i < 4; i++)
				cost += vBart(i)* vBar(i);
			
			pairs.at(i).cost = cost;
		}
		else
		{
			//pairs.at(i).vbar =  (mOBJ.m_pVertexList[edge.v1] + mOBJ.m_pVertexList[edge.v2])  * 0.5;
			pairs.at(i).vbar.push_back((mOBJ.m_pVertexList[edge.v1].x + mOBJ.m_pVertexList[edge.v2].x)  * 0.5);
			pairs.at(i).vbar.push_back((mOBJ.m_pVertexList[edge.v1].y + mOBJ.m_pVertexList[edge.v2].y)  * 0.5);
			pairs.at(i).vbar.push_back((mOBJ.m_pVertexList[edge.v1].z + mOBJ.m_pVertexList[edge.v2].z)  * 0.5);
			pairs.at(i).vbar.push_back(1);
			
			boost::numeric::ublas::vector <float> vBar (4);
			vBar(0) = pairs.at(i).vbar.at(0);
			vBar(0) = pairs.at(i).vbar.at(1);
			vBar(0) = pairs.at(i).vbar.at(2);
			vBar(0) = pairs.at(i).vbar.at(3);

		
			boost::numeric::ublas::vector <float> vBart (4);
			vBart  = prod(Q, vBar);
			//cout << vBar(0) << "  " << vBar(1) << "   "<<vBar(2) << "  " << vBar(3) << endl;
			float cost = 0;
			for (unsigned int i = 0; i < 4; i++)
				cost += vBart(i)* vBar(i);
			
			pairs.at(i).cost = cost;


		}
	
	}

	cout << "Q matrix initialization finished" << endl;
	int itr = 0;
	while (itr < MAX_ITR && mOBJ.m_nTriangles > 5 && mOBJ.m_nVertices > 5)
	{
		itr++;
	// sort pairs via cost
	sort(pairs.begin(),pairs.end(),compareCost);
#ifdef Debug
	//for (unsigned int i = 0; i < pairs.size();i++)
	//{
	//	cout << "v1:" << pairs.at(i).v1 << "  v2:" << pairs.at(i).v2 << endl;
	//}
	cout << "After sort according to cost = Q*vbar, itr = " << itr <<endl;
#endif
	// contract the top pair in pairs
	_pair current = pairs.at(0);
	mOBJ.m_pVertexList[current.v1] = Vec3f(current.vbar.at(0),  current.vbar.at(1),  current.vbar.at(2)) ;
	



	for (unsigned int i = 0; i < mOBJ.m_nTriangles; i++)
	{

		if (mOBJ.m_pTriangleList[i][0] == current.v2 || mOBJ.m_pTriangleList[i][1] == current.v2 || mOBJ.m_pTriangleList[i][2] == current.v2)
		{			
			// if [v1 v2] is an edge				
			if (mOBJ.m_pTriangleList[i][0] == current.v1 || mOBJ.m_pTriangleList[i][1] == current.v1 || mOBJ.m_pTriangleList[i][2] == current.v1)			
			{
				 for (unsigned int j = i; j < mOBJ.m_nTriangles -1 ; j++)
					 mOBJ.m_pTriangleList[j] = mOBJ.m_pTriangleList[j+1];
				
				mOBJ.m_nTriangles =  mOBJ.m_nTriangles - 1;
			}
			else
			{
				if (mOBJ.m_pTriangleList[i][0] == current.v2)
						mOBJ.m_pTriangleList[i][0] = current.v1;
				if (mOBJ.m_pTriangleList[i][1] == current.v2)
						mOBJ.m_pTriangleList[i][1] = current.v1;
				if (mOBJ.m_pTriangleList[i][2] == current.v2)
						mOBJ.m_pTriangleList[i][2] = current.v1;
			}
		}

	}

// %%%%%%%%%%%     update vertex list    %%%%%%%%%%

	// change vertex list and vertext that has an index > v2 in triangle list
		for (unsigned int i = current.v2 ; i < mOBJ.m_nVertices - 1 ; i++)
	{
		mOBJ.m_pVertexList[i] = mOBJ.m_pVertexList[i + 1];
	}
	mOBJ.m_nVertices--;

	// %%%%%%%%%        update triangle list: v2 first, v1 next      %%%%%%%%%%%%%%%%%%%%

	// update all index > v2			
			for (unsigned int i = current.v2 + 1; i < mOBJ.m_nVertices + 1 ; i++)
				{
					for (unsigned int j = 0; j < mOBJ.m_nTriangles; j++)
						{
							if (mOBJ.m_pTriangleList[j][0] == i)
								mOBJ.m_pTriangleList[j][0] = i - 1;
							if (mOBJ.m_pTriangleList[j][1] == i)
								mOBJ.m_pTriangleList[j][1] = i - 1;
							if (mOBJ.m_pTriangleList[j][2] == i)
								mOBJ.m_pTriangleList[j][2] = i - 1;
						}
				}			


	


	// %%%%%%%%%%if v1 does not exist in the current list, delete v1 and update all index > v1  %%%%%%%%%%%%%
	int cnt = 0;
	for (unsigned int j = 0; j < mOBJ.m_nTriangles ; j++)
		if (mOBJ.m_pTriangleList[j][0] == current.v1 || mOBJ.m_pTriangleList[j][1] == current.v1 || mOBJ.m_pTriangleList[j][2] == current.v1)
			{
				cnt ++;
			}

	if (cnt == 0)
	{
		for (unsigned int i = current.v1 ; i < mOBJ.m_nVertices - 1 ; i++)
		{
			mOBJ.m_pVertexList[i] = mOBJ.m_pVertexList[i + 1];
		}
		mOBJ.m_nVertices--;
		for (unsigned int i = current.v1 + 1; i < mOBJ.m_nVertices + 1 ; i++)
		{
				for (unsigned int j = 0; j < mOBJ.m_nTriangles; j++)
				{
						if (mOBJ.m_pTriangleList[j][0] == i)
															mOBJ.m_pTriangleList[j][0] = i - 1;
						if (mOBJ.m_pTriangleList[j][1] == i)
															mOBJ.m_pTriangleList[j][1] = i - 1;
						if (mOBJ.m_pTriangleList[j][2] == i)
															mOBJ.m_pTriangleList[j][2] = i - 1;
				}
		}

	}
	
	//for (unsigned int j = 0; j < mOBJ.m_nTriangles; j++)
	//{
	//		
	//	assert(mOBJ.m_pTriangleList[j][0] < mOBJ.m_nVertices);
	//	assert(mOBJ.m_pTriangleList[j][1] < mOBJ.m_nVertices);
	//	assert(mOBJ.m_pTriangleList[j][2] < mOBJ.m_nVertices);
	//}

	// %%%%%%%%%%%%%  update pairs   %%%%%%%%%%%%%%%%%%%%%
	// change pairs from v2 to v1
	pairs.erase(pairs.begin());
	for (unsigned int i = 0; i < pairs.size(); i++)
	{
		if (pairs.at(i).v2 == current.v2)
		{
			pairs.at(i).v2 == current.v1;
			if (pairs.at(i).v1 > pairs.at(i).v2)
			{
				unsigned int temp = pairs.at(i).v1;
				pairs.at(i).v1 = pairs.at(i).v2;
				pairs.at(i).v2 = temp;
			}
		}
		else if (pairs.at(i).v1 == current.v2)	
		{
			pairs.at(i).v1 == current.v1;
		}
				
			
	}
	for (unsigned int i = 0; i < pairs.size(); i++)
		assert(pairs.at(i).v1 != pairs.at(i).v2 );

	// %%%%%%%%%%%    update cost   %%%%%%%%%%%%%%%%%%
	for (unsigned int i = 0; i < pairs.size(); i++)
	{
		if (pairs.at(i).v1 == current.v1 ||  pairs.at(i).v2 == current.v1)	
		{
			_pair edge = pairs.at(i);
			boost::numeric::ublas::matrix <float> Q1 (4,4);
			Q1 = Qs[edge.v1];
		

			boost::numeric::ublas::matrix <float> Q2 (4,4);
			Q2 = Qs[edge.v2];
			boost::numeric::ublas::matrix <float> Q (4,4);
			Q = Q1 + Q2;
		
			boost::numeric::ublas::matrix <float> Q_inv (4,4);
			bool invertable = InvertMatrix(Q,Q_inv);
		
					if (invertable == true)
		{
			
			pairs.at(i).vbar.push_back(Q_inv(0,3) / Q_inv(3,3));
			pairs.at(i).vbar.push_back(Q_inv(1,3) / Q_inv(3,3));
			pairs.at(i).vbar.push_back(Q_inv(2,3) / Q_inv(3,3));
			pairs.at(i).vbar.push_back(1);
			
			
			boost::numeric::ublas::vector <float> vBar (4);
			vBar(0) = pairs.at(i).vbar.at(0);
			vBar(0) = pairs.at(i).vbar.at(1);
			vBar(0) = pairs.at(i).vbar.at(2);
			vBar(0) = pairs.at(i).vbar.at(3);

			boost::numeric::ublas::vector <float> vBart (4);
			vBart  = prod(Q, vBar);
			//cout << vBar(0) << "  " << vBar(1) << "   "<<vBar(2) << "  " << vBar(3) << endl;
			float cost = 0;
			for (unsigned int i = 0; i < 4; i++)
				cost += vBart(i)* vBar(i);
			
			pairs.at(i).cost = cost;
		}
		else
		{
			//pairs.at(i).vbar =  (mOBJ.m_pVertexList[edge.v1] + mOBJ.m_pVertexList[edge.v2])  * 0.5;
			pairs.at(i).vbar.push_back((mOBJ.m_pVertexList[edge.v1].x + mOBJ.m_pVertexList[edge.v2].x)  * 0.5);
			pairs.at(i).vbar.push_back((mOBJ.m_pVertexList[edge.v1].y + mOBJ.m_pVertexList[edge.v2].y)  * 0.5);
			pairs.at(i).vbar.push_back((mOBJ.m_pVertexList[edge.v1].z + mOBJ.m_pVertexList[edge.v2].z)  * 0.5);
			pairs.at(i).vbar.push_back(1);
			
			boost::numeric::ublas::vector <float> vBar (4);
			vBar(0) = pairs.at(i).vbar.at(0);
			vBar(0) = pairs.at(i).vbar.at(1);
			vBar(0) = pairs.at(i).vbar.at(2);
			vBar(0) = pairs.at(i).vbar.at(3);

		
			boost::numeric::ublas::vector <float> vBart (4);
			vBart  = prod(Q, vBar);
			//cout << vBar(0) << "  " << vBar(1) << "   "<<vBar(2) << "  " << vBar(3) << endl;
			float cost = 0;
			for (unsigned int i = 0; i < 4; i++)
				cost += vBart(i)* vBar(i);
			
			pairs.at(i).cost = cost;


		}

		}
	}
				
	}    // end while (itr < MAXITERATION)
	
	for (unsigned int i = 0; i < mOBJ.m_nTriangles; i++)
	{
		Vec3f v0 = mOBJ.m_pVertexList[mOBJ.m_pTriangleList[i][0]];
		Vec3f v1 = mOBJ.m_pVertexList[mOBJ.m_pTriangleList[i][1]];
		Vec3f v2 = mOBJ.m_pVertexList[mOBJ.m_pTriangleList[i][2]];
		if ((v1-v2).L2Norm_Sqr() > MAX_POINT_DIS || (v0-v2).L2Norm_Sqr() > MAX_POINT_DIS  || (v1-v0).L2Norm_Sqr() > MAX_POINT_DIS )
		{
			for (unsigned int j = i; j < mOBJ.m_nTriangles -1 ;j++)
					 mOBJ.m_pTriangleList[j] = mOBJ.m_pTriangleList[j+1];
			mOBJ.m_nTriangles =  mOBJ.m_nTriangles - 1;
		}
		if (abs (mOBJ.m_pTriangleList[i][0] - mOBJ.m_pTriangleList[i][1]) > 100 || abs (mOBJ.m_pTriangleList[i][0] - mOBJ.m_pTriangleList[i][2]) > 100 
			|| abs (mOBJ.m_pTriangleList[i][1] - mOBJ.m_pTriangleList[i][2]) > 100 )
			{
				for (unsigned int j = i; j < mOBJ.m_nTriangles -1 ;j++)
					 mOBJ.m_pTriangleList[j] = mOBJ.m_pTriangleList[j+1];
			mOBJ.m_nTriangles =  mOBJ.m_nTriangles - 1;
		}

		Vec3f aa = v0 - v1;
		Vec3f bb = v1 - v2;
		Vec3f cc = v0 - v2;
		aa.Normalize();
		bb.Normalize();
		cc.Normalize();
		if ( (aa.x * bb.x +  aa.y * bb.y + aa.z * bb.z) > 0.9  || (aa.x * cc.x +  aa.y * cc.y + cc.z * bb.z)> 0.9 || (cc.x * bb.x +  cc.y * bb.y + cc.z * bb.z)> 0.9)

		{
			for (unsigned int j = i; j < mOBJ.m_nTriangles -1 ;j++)
					 mOBJ.m_pTriangleList[j] = mOBJ.m_pTriangleList[j+1];
			mOBJ.m_nTriangles =  mOBJ.m_nTriangles - 1;
		}

	}


	cout << "the final number of vertices is " << mOBJ.m_nVertices << endl; 
	cout << "the final number of triangles is " << mOBJ.m_nTriangles << endl; 
	
	SimpleOBJ::CSimpleObject test;

	bool m_loadSuccess =  test.LoadFromObj(fileName);
	float totaldis = 0;
	for (unsigned int i = 0; i < test.m_nVertices; i++)
	{
		float MinDis = 100000;
		for (unsigned int j = 0; j < mOBJ.m_nVertices; j++)
			if (      (test.m_pVertexList[i] - mOBJ.m_pVertexList[j]).L2Norm_Sqr() < MinDis )
			{
				MinDis = (test.m_pVertexList[i] - mOBJ.m_pVertexList[j]).L2Norm_Sqr();
			}
		totaldis += MinDis;
	}

	totaldis /= mOBJ.m_nVertices;
	 totaldis /= test.m_nVertices;
	 cout << "the evaluation score is " <<  totaldis<< endl; 
	mOBJ.SaveToObj("simp.obj");
	return 0;
}

