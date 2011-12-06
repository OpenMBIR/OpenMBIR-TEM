#ifndef _WrappedArray_hpp_
#define _WrappedArray_hpp_


#include "MXA/Common/MXASetGetMacros.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/allocate.h"

template< typename T>
class Wrapped3DArray
{
  public:
      typedef WrappedArray<T>                      Self;
      typedef boost::shared_ptr<Self >        Pointer;
      typedef boost::shared_ptr<const Self >  ConstPointer;
      typedef boost::weak_ptr<WrappedArray<T> > WeakPointer;
      typedef boost::weak_ptr<WrappedArray<T> > ConstWeakPointer;
      static Pointer NullPointer(void)
      {
        return Pointer(static_cast<WrappedArray<T>*>(0));
      }


      static Pointer New(size_t Dim0=0, size_t Dim1=0, size_t Dim2=0)
      {
        Pointer sharedPtr (new Wrapped3DArray<T>(Dim0, Dim1, Dim2));
        return sharedPtr;
      }

    virtual ~Wrapped3DArray()
    {
      for (size_t i = 0; i < m_Dim1; ++i)
      {
        for (size_t j = 0; j < m_Dim2; ++j)
        {
           delete [] ptr[i][j];
        }
        delete [] ptr[i];
      }
      delete [] ptr;
    }

    T*** getPointer() { return ptr; }

  protected:
    Wrapped3DArray(size_t Dim0=0, size_t Dim1=0, size_t Dim2=0) :
      m_Dim0(Dim0),
      m_Dim1(Dim1),
      m_Dim2(Dim2)
    {
//      if (Dim0 * Dim1 * Dim2 == 0) { ptr = NULL; return; }
//      ptr = new T*[Dim0];
//      for(size_t i = 0; i < Dim1; ++i)
//      {
//        ptr[i] = new T[Dim1];
//        for (size_t j = 0; j < Dim2; ++j)
//        {
//          ptr[i][j] = new T[Dim2];
//        }
//      }
    }



  private:
    T*** ptr;
    size_t m_Dim0;
    size_t m_Dim1;
    size_t m_Dim2;

    Wrapped3DArray<(const Wrapped3DArray&); // Copy Constructor Not Implemented
    void operator=(const Wrapped3DArray&); // Operator '=' Not Implemented
};




#endif /* _WrappedArray_hpp_  */
