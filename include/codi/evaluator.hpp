/**
 * Evaluator class (draft by J. Fuhrmann <juergen.fuhrmann@wias-berlin.de>)
 * This has been inspired by JuliaDiff/DiffResults.jl. The main purpose is to provide a convenient way
 * to evaluate function and jacobian at once, as it is often used in Newton solvers.
 **/


#include <exception>
#include <vector>

namespace codi
{

  /**
   * @brief Evaluator of value and jacobian of vector function with origin dimension as template parameter
   *
   *
   * It takes as template argument a function f(x,y) mapping an origin vector x of dimension nx to
   * an image vector y of dimension ny and allows to evaluate function and derivative.
   *
   * The vector data type is variable, by default, std::vector is used.
   */
  template <
    class FUNC, // Function  R^nx -> R^ny
    int NX,     // Fixed  origin dimension
    template<class T, class ALLOC> class VECTOR=std::vector
    >
  class FixedDimensionEvaluator
  {
    typedef RealForwardVec<NX> codivec;  // CODI vector type

  private:
    FUNC func;      // Function to be evaluated
    const int nx;   // Origin dimension
    const int ny;   // Image dimension
    VECTOR<codivec,std::allocator<codivec>> x; // Vector of unknown values
    VECTOR<codivec,std::allocator<codivec>> y; // Vector of result values
    
  public:
    FixedDimensionEvaluator(int nx, int ny, FUNC& func):
      func(func),
      nx(nx),
      ny(ny),
      x(nx),
      y(ny)
    {
      if (nx!=NX)
        throw std::domain_error("codi::FixedDimensionEvaluator: template parameter NX must coincide with origin dimension nx.");
    };

    /**
     * @brief Call function for vector x
     *
     * Note: we might check for matching dimensions
     */
    template <class V>
    void call(const V &x)
    {
      // Initialize x and gradient
      for (int i=0;i<nx;i++)
      {
        this->x[i]=x[i];
        this->x[i].gradient()[i]=1.0;
      }
      // Call function 
      func(this->x,this->y);
    }

    
    /**
     * @brief Retrieve function value f_i(x)
     *
     * Note: we might check for matching dimensions
     */
    double result(int i)
    {
      return y[i].getValue();
    }

    /**
     * @brief Retrieve jacobian value at d f_i(x)/d x_j
     *
     * Note: we might check for matching dimensions
     */
    double jacobian(int i, int j)
    {
      return y[i].getGradient()[j];
    }
  };
  
  
  
  /**
   * @brief Evaluator of value and jacobian of vector function with arbitrary origin dimension
   *
   * It takes as template argument a function f(x,y) mapping an origin vector x of dimension nx to
   * an image vector y of dimension ny and allows to evaluate function and derivative.
   *
   * The vector data type is variable, by default, std::vector is used.
   */
  template <
    class FUNC,  // Function  R^nx -> R^ny
    template <class T, class ALLOC> class VECTOR=std::vector // vector type
    >
  class Evaluator
  {
    const int nx;

    // Work around  fixed origin dimension: create a union of FixedDimensionEvaluators
    // with different template dimensions
    union
    {
#define CODE(i) FixedDimensionEvaluator<FUNC,i, VECTOR> * eval##i;
#include "codi/evaluator-macros.h"
#undef CODE
    } eval;
    
  public:
    Evaluator(int nx, int ny, FUNC& func):
      nx(nx)
    {

      // call constructor for origin dimension nx
      switch(nx)
      {
#define CODE(i) case i: eval.eval##i=new FixedDimensionEvaluator<FUNC,i,VECTOR>(nx,ny,func); break;
#include "codi/evaluator-macros.h"
#undef CODE
      default:
        std::stringstream msg;
        msg << "codi::Evaluator: unable to handle dimensions larger than " << CODI_MAXDIM;
        throw std::domain_error(msg.str());
      }
      
    }
    
    
    /**
     * @brief Call function for vector x
     *
     * The call is delegated to the proper FixedDimensionEvaluator
     */
    template <class V>
    void call(const V& x)
    {
      // Delegate call
      switch(nx)
      {
#define CODE(i) case i: eval.eval##i->call(x); break;
#include "codi/evaluator-macros.h"
#undef CODE
      default: break;
      }
    }
    
    /**
     * @brief Retrieve function value f_i(x)
     *
     * Note: we might check for matching dimensions
     */
    double result(int iy)
    {
      // Delegate call 
      switch(nx)
      {
#define CODE(i) case i: return eval.eval##i->result(iy); 
#include "codi/evaluator-macros.h"
#undef CODE
      default: return 0;
      }
    }
    
    /**
     * @brief Retrieve jacobian value at d f_i(x)/d x_j
     *
     * Note: we might check for matching dimensions
     */
    double jacobian(int iy, int ix)
    {
      // Delegate call 
      switch(nx)
      {
#define CODE(i) case i: return eval.eval##i->jacobian(iy,ix); 
#include "codi/evaluator-macros.h"
#undef CODE
      default: return 0;
      }
    }
    
    ~Evaluator()
    {
      // Destroy the right FixedDimensionEvaluator
      switch(nx)
      {
#define CODE(i) case i: delete eval.eval##i; break;
#include "codi/evaluator-macros.h"
#undef CODE
      default: break;
      }
    }
  };
  

#if __cplusplus >= 201703L
  // C++ 17 template argument deduction allows to avoid template
  // parameter in constructor of Evaluator
  template<
    class FUNC,
    template <class T, class ALLOC> class VECTOR=std::vector
    >
  Evaluator(int x, int y, FUNC&func)  -> Evaluator<FUNC,VECTOR>;
#endif
  
}
