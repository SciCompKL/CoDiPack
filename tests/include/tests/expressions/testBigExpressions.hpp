#include "../../testInterface.hpp"

struct TestBigExpressions : public TestInterface {
  public:
    NAME("BigExpressions")
    IN(12)
    OUT(1)
    POINTS(1) =  // clang-format off
    {
      {1.25, 2.5, 3.25, 4.5, 5.75, 6.25, 7.5, 8.5, 9.25, 10.25, 11.75, 12.5}
    
    };  // clang-format on
    
    template<typename Number>
    static void func(Number* x, Number* y) {
      Number* val_coord_Edge_CG = &x[0];
      Number* val_coord_FaceElem_CG = &x[3];
      Number* val_coord_Elem_CG = &x[6];
      Number* val_coord_Point = &x[9];
    
      auto vec_a0 =  val_coord_Edge_CG[0]-val_coord_Point[0];
      auto vec_a1 =  val_coord_Edge_CG[1]-val_coord_Point[1];
      auto vec_a2 =  val_coord_Edge_CG[2]-val_coord_Point[2];
    
      auto vec_b0 =  val_coord_FaceElem_CG[0]-val_coord_Point[0];
      auto vec_b1 =  val_coord_FaceElem_CG[1]-val_coord_Point[1];
      auto vec_b2 =  val_coord_FaceElem_CG[2]-val_coord_Point[2];
    
      auto vec_c0 =  val_coord_Elem_CG[0]-val_coord_Point[0];
      auto vec_c1 =  val_coord_Elem_CG[1]-val_coord_Point[1];
      auto vec_c2 =  val_coord_Elem_CG[2]-val_coord_Point[2];
    
      auto vec_d0 = vec_a1*vec_b2-vec_a2*vec_b1;
      auto vec_d1 = -(vec_a0*vec_b2-vec_a2*vec_b0);
      auto vec_d2 = vec_a0*vec_b1-vec_a1*vec_b0;
    //  Number vec_a0 =  val_coord_Edge_CG[0]-val_coord_Point[0];
    //  Number vec_a1 =  val_coord_Edge_CG[1]-val_coord_Point[1];
    //  Number vec_a2 =  val_coord_Edge_CG[2]-val_coord_Point[2];
    //
    //  Number vec_b0 =  val_coord_FaceElem_CG[0]-val_coord_Point[0];
    //  Number vec_b1 =  val_coord_FaceElem_CG[1]-val_coord_Point[1];
    //  Number vec_b2 =  val_coord_FaceElem_CG[2]-val_coord_Point[2];
    //
    //  Number vec_c0 =  val_coord_Elem_CG[0]-val_coord_Point[0];
    //  Number vec_c1 =  val_coord_Elem_CG[1]-val_coord_Point[1];
    //  Number vec_c2 =  val_coord_Elem_CG[2]-val_coord_Point[2];
    //
    //  Number vec_d0 = vec_a1*vec_b2-vec_a2*vec_b1;
    //  Number vec_d1 = -(vec_a0*vec_b2-vec_a2*vec_b0);
    //  Number vec_d2 = vec_a0*vec_b1-vec_a1*vec_b0;
    
      y[0] = fabs(vec_c0*vec_d0 + vec_c1*vec_d1 + vec_c2*vec_d2)/6.0;
    
    }
};
