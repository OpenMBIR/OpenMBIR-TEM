#ifndef _VoxelUpdateList_H_
#define _VoxelUpdateList_H_

#include <ostream>

#include <boost/shared_array.hpp>


#include "MXA/Common/MXASetGetMacros.h"


#include "MBIRLib/MBIRLib.h"

/**
 * @brief The VoxelUpdateList class holds the x & z index values for a voxel update list. The class uses the boost::shared_array
 * to automatically clean up the array when all references to the class go out of scope. The data is laid out in memory
 * such that the x & z values are side by side (interleaved)
 */
class MBIRLib_EXPORT VoxelUpdateList
{
  public:
    MXA_SHARED_POINTERS(VoxelUpdateList)
    MXA_TYPE_MACRO(VoxelUpdateList)

    typedef boost::shared_array<int32_t> ArrayType;


    /**
     * @brief New
     * @param numElements
     * @return
     */
    static Pointer New(int32_t numElements);

    /**
     * @brief New Creates a new VoxelUpdateList that uses an externally created array as its data source
     * @param numElements
     * @param array
     * @return
     */
    static Pointer New(int32_t numElements, ArrayType array);

    /**
     * @brief GenRandList
     * @param InList
     * @return
     */
    static Pointer GenRandList(Pointer InList, bool print = false);

    /**
     * @brief GenRegularList
     * @param jCount
     * @return
     */
    static Pointer GenRegularList(uint16_t jCount, uint16_t kCount);

    /**
     * @brief ~VoxelUpdateList
     */
    virtual ~VoxelUpdateList();

    /**
     * @brief numElements
     * @return
     */
    int32_t numElements();

    /**
     * @brief setNumElements
     */
    void setNumElements();

    /**
     * @brief xIdx
     * @param index
     * @return
     */
    int32_t xIdx(int32_t index);
    void setX(int32_t index, int32_t x);

    /**
     * @brief zIdx
     * @param index
     * @return
     */
    int32_t zIdx(int32_t index);
    void setZ(int32_t index, int32_t z);

    /**
     * @brief setPair
     * @param index
     * @param x
     * @param z
     */
    void setPair(int32_t index, int32_t x, int32_t z);

    /**
     * @brief resize
     * @param numElements
     */
    void resize(int32_t numElements);

    /**
     * @brief getArray
     * @return
     */
    ArrayType getArray();

    /**
     * @brief subList
     * @param numElements
     * @param start
     * @return
     */
    Pointer subList(int32_t numElements, int32_t start);

    /**
     * @brief printList
     * @param out
     */
    void printList(std::ostream& out);
    /**
     * @brief printMaxList
     * @param out
     */
    void printMaxList(std::ostream& out);
    /**
     * @brief printMinList
     * @param out
     */
    void printMinList(std::ostream& out);



  protected:
    VoxelUpdateList(int32_t numElements);

  private:
    int32_t  m_NumElements;
    ArrayType m_Array;

    VoxelUpdateList(const VoxelUpdateList&); // Copy Constructor Not Implemented
    void operator=(const VoxelUpdateList&); // Operator '=' Not Implemented

};



#endif /* _VoxelUpdateList_H_ */
