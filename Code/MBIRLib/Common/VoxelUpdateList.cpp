#include "MBIRLib/Common/VoxelUpdateList.h"

#include <limits>

//-- Boost Headers for Random Numbers
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include "MBIRLib/Common/EIMTime.h"
#include "MBIRLib/Reconstruction/ReconstructionStructures.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdateList::VoxelUpdateList(int32_t numElements) :
  m_NumElements(numElements)
{
  m_Array = ArrayType(new int32_t[m_NumElements * 2]);
  ::memset(m_Array.get(), 0, m_NumElements * 2); // Initialize everthing to 0
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdateList::Pointer VoxelUpdateList::New(int32_t numElements)
{
  VoxelUpdateList::Pointer ptr(new VoxelUpdateList(numElements) );
  return ptr;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdateList::Pointer VoxelUpdateList::New(int32_t numElements, ArrayType array)
{
  Pointer ptr = VoxelUpdateList::New(numElements);
  ptr->m_Array = array;
  return ptr;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdateList::~VoxelUpdateList()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t VoxelUpdateList::numElements()
{
  return m_NumElements;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t VoxelUpdateList::xIdx(int32_t index)
{
  BOOST_ASSERT(index < (m_NumElements) );
  return m_Array[index * 2];
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void VoxelUpdateList::setX(int32_t index, int32_t x)
{
  BOOST_ASSERT(index < (m_NumElements) );
  m_Array[index * 2] = x;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t VoxelUpdateList::zIdx(int32_t index)
{
  BOOST_ASSERT(index < (m_NumElements) );
  return m_Array[index * 2 + 1];
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void VoxelUpdateList::setZ(int32_t index, int32_t z)
{
  BOOST_ASSERT(index < (m_NumElements) );

  m_Array[index * 2 + 1]= z;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void VoxelUpdateList::setPair(int32_t index, int32_t x, int32_t z)
{
  BOOST_ASSERT(index < (m_NumElements) );
  m_Array[index * 2] = x;
  m_Array[index * 2 + 1] = z;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void VoxelUpdateList::resize(int32_t numElements)
{
  if (numElements == m_NumElements) { return; }
  if (numElements <= 0)
  {
    m_NumElements = 0;
    m_Array.reset(); // This will destroy the internal array
    return;
  }
  ArrayType data = ArrayType(new int32_t[numElements * 2]);
  ::memset(data.get(), 0, numElements * 2); // Initialize everthing to 0

  // Copy the data from the old array into the new array
  if (numElements > m_NumElements) // The new array is larger than the older array
  {
    ::memcpy(data.get(), m_Array.get(), 2 * m_NumElements * sizeof(int32_t) );
  }
  else // The new array is SMALLER than the older array
  {
    ::memcpy(data.get(), m_Array.get(), 2 * numElements * sizeof(int32_t) );
  }
  m_Array = data;
  m_NumElements = numElements;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdateList::Pointer VoxelUpdateList::subList(int32_t numElements, int32_t start)
{
  BOOST_ASSERT(start + numElements <= m_NumElements);
  Pointer ptr = VoxelUpdateList::New(numElements);
  // Copy our data into the new array
  ::memcpy(ptr->getArray().get(), &(m_Array[start * 2]), numElements * 2 * sizeof(int32_t) );
  return ptr;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdateList::ArrayType VoxelUpdateList::getArray()
{
  return m_Array;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdateList::Pointer VoxelUpdateList::GenRandList(VoxelUpdateList::Pointer InList)
{
  // Create a Copy
  VoxelUpdateList::Pointer OpList = VoxelUpdateList::New(InList->numElements());

  const uint32_t rangeMin = 0;
  const uint32_t rangeMax = std::numeric_limits<uint32_t>::max();
  typedef boost::uniform_int<uint32_t> NumberDistribution;
  typedef boost::mt19937 RandomNumberGenerator;
  typedef boost::variate_generator<RandomNumberGenerator&,NumberDistribution> Generator;

  NumberDistribution distribution(rangeMin, rangeMax);
  RandomNumberGenerator generator;
  Generator numberGenerator(generator, distribution);
  boost::uint32_t arg = static_cast<boost::uint32_t>(EIMTOMO_getMilliSeconds());
  generator.seed(arg); // seed with the current time


  size_t dims[3] =
  { InList->numElements(), 0, 0};
  Int32ArrayType::Pointer Counter = Int32ArrayType::New(dims, "Counter");

  for (int32_t j_new = 0; j_new < InList->numElements(); j_new++)
  {
    Counter->d[j_new] = j_new;
  }

  int32_t ArraySize = InList->numElements();
  ArrayType OpListPtr = OpList->getArray();
  ArrayType InptListPtr = InList->getArray();

  for(int32_t test_iter = 0;test_iter < InList->numElements();test_iter++)
  {
    int32_t Index = numberGenerator() % InList->numElements();
    OpListPtr[test_iter] = InptListPtr[Counter->d[Index]];
    Counter->d[Index] = Counter->d[InList->numElements() - 1];
    ArraySize--;
  }

#ifdef DEBUG
  std::cout<<"Input"<<std::endl;
  InList->printMaxList(std::cout);

  std::cout<<"Output"<<std::endl;
  OpList->printMaxList(std::cout);
#endif //Debug

  // Return the Randomized List
  return OpList;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
VoxelUpdateList::Pointer VoxelUpdateList::GenRegularList(uint16_t jCount, uint16_t kCount)
{
  Pointer InpList = New(jCount * kCount);

  int32_t iter = 0;
  for (int16_t j = 0; j < jCount; j++)
  {
    for (int16_t k = 0; k < kCount; k++)
    {
      InpList->setPair(iter, k, j);
      iter++;
    }
  }
  return InpList;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void VoxelUpdateList::printList(std::ostream &out)
{
  out<<"Printing List :"<<std::endl;
  for(int32_t i=0; i < m_NumElements; i++)
    out<<"("<<m_Array[i*2 + 1]<<","<<m_Array[i*2]<<")"<<std::endl;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void VoxelUpdateList::printMaxList(std::ostream &out)
{
  int32_t max_z = std::numeric_limits<int32_t>::min();
  int32_t max_x = std::numeric_limits<int32_t>::min();
  for(int32_t i = 0; i < m_NumElements; i++)
  {
    if(m_Array[i*2 + 1] > max_z)
      max_z = m_Array[i*2 + 1];

    if(m_Array[i*2] > max_x)
      max_x = m_Array[i*2];
  }
  out<<"Number of elements ="<<m_NumElements<<std::endl;
  out<<"("<<max_z<<","<<max_x<<")"<<std::endl;
}
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void VoxelUpdateList::printMinList(std::ostream &out)
{
  int32_t min_z = std::numeric_limits<int32_t>::max(); //pow(2, 31);
  int32_t min_x = std::numeric_limits<int32_t>::max();
  for(int32_t i = 0; i < m_NumElements; i++)
  {
    if(m_Array[i*2 + 1] < min_z)
      min_z = m_Array[i*2 + 1];

    if(m_Array[i*2] < min_x)
      min_x = m_Array[i*2];
  }
  out<<"Number of elements ="<<m_NumElements<<std::endl;
  out<<"("<<min_z<<","<<min_x<<")"<<std::endl;
}

