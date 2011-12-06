/*
 * CalculateAMatrix.h
 *
 *  Created on: Dec 5, 2011
 *      Author: mjackson
 */

#ifndef CALCULATEAMATRIX_H_
#define CALCULATEAMATRIX_H_

#include "MXA/Common/MXASetGetMacros.h"

#include "TomoEngine/TomoEngine.h"
#include "TomoEngine/Common/AbstractFilter.h"
#include "TomoEngine/SOC/SOCStructures.h"

/*
 *
 */
class CalculateAMatrixColumn : public AbstractFilter
{
  public:

    MXA_SHARED_POINTERS(CalculateAMatrixColumn);
    MXA_STATIC_NEW_MACRO(CalculateAMatrixColumn);
    MXA_TYPE_MACRO_SUPER(CalculateAMatrixColumn, AbstractFilter);

    virtual ~CalculateAMatrixColumn();

    MXA_INSTANCE_PROPERTY(TomoInputs*, Inputs);
    MXA_INSTANCE_PROPERTY(Sinogram*, Sinogram);
    MXA_INSTANCE_PROPERTY(Geometry*, Geometry);

    MXA_INSTANCE_PROPERTY_OLD(DATA_TYPE*, Cosine, cosine);
    MXA_INSTANCE_PROPERTY_OLD(DATA_TYPE*, Sine, sine);

    MXA_INSTANCE_PROPERTY_OLD(DATA_TYPE, BeamWidth, BEAM_WIDTH);
    MXA_INSTANCE_PROPERTY_OLD(uint16_t, row, row);
    MXA_INSTANCE_PROPERTY_OLD(uint16_t, col, col);
    MXA_INSTANCE_PROPERTY_OLD(uint16_t, Slice, slice);
    MXA_INSTANCE_PROPERTY_OLD(DATA_TYPE**, VoxelProfile, VoxelProfile);
    MXA_INSTANCE_PROPERTY_OLD(DATA_TYPE*, D1, d1);
    MXA_INSTANCE_PROPERTY_OLD(DATA_TYPE*, D2, d2);
    MXA_INSTANCE_PROPERTY_OLD(AMatrixCol*, AMatrixCol, Ai)


    virtual void execute();

  protected:
    CalculateAMatrixColumn();
};

#endif /* CALCULATEAMATRIX_H_ */
