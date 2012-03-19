#error
uint16_t subIterations = 1;

if (updateType == RegularRandomOrderUpdate)
{
  std::cout << indent << "Regular Random Order update of Voxels" << std::endl;
}
else if (updateType == HomogeniousUpdate)
{
  std::cout << indent << "Homogenous update of voxels" << std::endl;
}
else if (updateType == NonHomogeniousUpdate)
{
  std::cout << indent << "Non Homogenous update of voxels" << std::endl;
  subIterations = NUM_NON_HOMOGENOUS_ITER;
}
else
{
  std::cout << indent << "Unknown Voxel Update Type. Returning Now" << std::endl;
  return;
}

DATA_TYPE NH_Threshold = 0.0;

for (uint16_t NH_Iter = 0; NH_Iter < subIterations; ++NH_Iter)
{
  if (updateType == NonHomogeniousUpdate)
  {
    //Compute VSC and create a map of pixels that are above the threshold value
    ComputeVSC();
    NH_Threshold = SetNonHomThreshold();
    //Use  FiltMagUpdateMap  to find MagnitudeUpdateMask
    //std::cout << "Completed Calculation of filtered magnitude" << std::endl;
    //Calculate the threshold for the top ? % of voxel updates
  }

  //printf("Iter %d\n",Iter);
#ifdef ROI
  AverageUpdate = 0;
  AverageMagnitudeOfRecon = 0;
#endif

#ifdef RANDOM_ORDER_UPDATES
  ArraySize = m_Geometry->N_x * m_Geometry->N_z;
  for (j_new = 0; j_new < ArraySize; j_new++)
  {
    Counter->d[j_new] = j_new;
  }
  for (j = 0; j < m_Geometry->N_z; j++)
  {
    for (k = 0; k < m_Geometry->N_x; k++)
    {
      if (updateType == NonHomogeniousUpdate)
      {
        // WHERE DOES MagUpdateMap get initialized at? These if statements will not cover all the possibilities
        if(MagUpdateMap->d[j][k] > NH_Threshold)
        {
          MagUpdateMask->d[j][j] = 1;
          MagUpdateMap->d[j][k] = 0;
        }
        else
        {
          MagUpdateMask->d[j][j] = 0;
        }
      }
      else if (updateType == HomogeniousUpdate)
      {
        MagUpdateMap->d[j][k] = 0;
      }

      VisitCount->d[j][k] = 0;
    }
  }
#endif

  START_TIMER;
  for (j = 0; j < m_Geometry->N_z; j++) //Row index
  {
    for (k = 0; k < m_Geometry->N_x; k++) //Column index
    {
#ifdef RANDOM_ORDER_UPDATES
      //RandomNumber=init_genrand(Iter);
      Index = (genrand_int31(RandomNumber)) % ArraySize;
      k_new = Counter->d[Index] % m_Geometry->N_x;
      j_new = Counter->d[Index] / m_Geometry->N_x;
      //memmove(Counter+Index,Counter+Index+1,sizeof(int32_t)*(ArraySize - Index-1));
      //TODO: Instead just swap the value in Index with the one in ArraySize
      Counter->d[Index] = Counter->d[ArraySize - 1];
      VisitCount->d[j_new][k_new] = 1;
      ArraySize--;
#else
      j_new=j;
      k_new=k;
#endif //Random order updates
      TempMemBlock = TempCol[j_new][k_new]; //Remove this

      int shouldInitNeighborhood = 0;
      if(TempMemBlock->count > 0) { ++shouldInitNeighborhood;}
      if (updateType == NonHomogeniousUpdate && shouldInitNeighborhood == 1 && MagUpdateMask->d[j_new][k_new] == 1) { ++shouldInitNeighborhood;}

      if(shouldInitNeighborhood > 0 )
      //After this should ideally call UpdateVoxelLine(j_new,k_new) ie put everything in this "if" inside a method called UpdateVoxelLine
      {
        for (i = 0; i < m_Geometry->N_y; i++) //slice index
        {
          //Neighborhood of (i,j,k) should be initialized to zeros each time
          for (int32_t p = 0; p <= 2; p++)
          {
            for (int32_t q = 0; q <= 2; q++)
            {
              for (r = 0; r <= 2; r++)
              {
                NEIGHBORHOOD[p][q][r] = 0.0;
                BOUNDARYFLAG[p][q][r] = 0;
              }
            }
          }
#ifdef CIRCULAR_BOUNDARY_CONDITION
          for(p = -1; p <=1; p++)
          for(q = -1; q <= 1; q++)
          for(r = -1; r <= 1;r++)
          {
            tempindex_x = mod(r+k_new,m_Geometry->N_x);
            tempindex_y =mod(p+i,m_Geometry->N_y);
            tempindex_z = mod(q+j_new,m_Geometry->N_z);
            NEIGHBORHOOD[p+1][q+1][r+1] = m_Geometry->Object->d[tempindex_z][tempindex_x][tempindex_y];
            BOUNDARYFLAG[p+1][q+1][r+1]=1;
          }
#else
          //For a given (i,j,k) store its 26 point neighborhood
          for (int32_t p = -1; p <= 1; p++)
          {
            for (int32_t q = -1; q <= 1; q++)
            {
              for (r = -1; r <= 1; r++)
              {
                if(i + p >= 0 && i + p < m_Geometry->N_y)
                {
                  if(j_new + q >= 0 && j_new + q < m_Geometry->N_z)
                  {
                    if(k_new + r >= 0 && k_new + r < m_Geometry->N_x)
                    {
                      NEIGHBORHOOD[p + 1][q + 1][r + 1] = m_Geometry->Object->d[q + j_new][r + k_new][p + i];
                      BOUNDARYFLAG[p + 1][q + 1][r + 1] = 1;
                    }
                    else
                    {
                      BOUNDARYFLAG[p + 1][q + 1][r + 1] = 0;
                    }
                  }
                }
              }
            }
          }
#endif//circular boundary condition check
          NEIGHBORHOOD[1][1][1] = 0.0;

#ifndef NDEBUG
          if(i == 0 && j == 31 && k == 31)
          {
            printf("***************************\n");
            printf("Geom %lf\n", m_Geometry->Object->d[i][31][31]);
            for (int p = 0; p <= 2; p++)
            {
              for (int q = 0; q <= 2; q++)
              {
                for (r = 0; r <= 2; r++)
                {
                  printf("%lf\n", NEIGHBORHOOD[p][q][r]);
                }
              }
            }
          }
#endif
          //Compute theta1 and theta2
          V = m_Geometry->Object->d[j_new][k_new][i];//Store the present value of the voxel
          THETA1 = 0.0;
          THETA2 = 0.0;

          //TempCol = CE_CalculateAMatrixColumn(j, k, i, Sinogram, Geometry, VoxelProfile);
          for (uint32_t q = 0; q < TempMemBlock->count; q++)
          {
            uint16_t i_theta = floor(static_cast<float>(TempMemBlock->index[q] / (m_Sinogram->N_r)));
            uint16_t i_r = (TempMemBlock->index[q] % (m_Sinogram->N_r));
            VoxelLineAccessCounter = 0;
            for (uint32_t i_t = VoxelLineResponse[i].index[0]; i_t < VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count; i_t++)
            {
              THETA2 += ((NuisanceParams.I_0->d[i_theta] * NuisanceParams.I_0->d[i_theta])
                  * (VoxelLineResponse[i].values[VoxelLineAccessCounter] * VoxelLineResponse[i].values[VoxelLineAccessCounter])
                  * (TempMemBlock->values[q]) * (TempMemBlock->values[q]) * Weight->d[i_theta][i_r][i_t]);
              THETA1 += NuisanceParams.I_0->d[i_theta] * ErrorSino->d[i_theta][i_r][i_t] * (TempMemBlock->values[q])
              * (VoxelLineResponse[i].values[VoxelLineAccessCounter]) * Weight->d[i_theta][i_r][i_t];
              VoxelLineAccessCounter++;
            }
          }
          THETA1 *= -1;
          minMax(&low, &high);

#ifdef DEBUG
          if(i == 0 && j == 31 && k == 31) printf("(%lf,%lf,%lf) \n", low, high, V - (THETA1 / THETA2));
#endif

          //Solve the 1-D optimization problem
          //printf("V before updating %lf",V);
#ifndef SURROGATE_FUNCTION
          //TODO : What if theta1 = 0 ? Then this will give error
          DerivOfCostFunc docf(BOUNDARYFLAG, NEIGHBORHOOD, FILTER, V, THETA1, THETA2, SIGMA_X_P, MRF_P);

          UpdatedVoxelValue = (DATA_TYPE)solve<DerivOfCostFunc>(&docf, (double)low, (double)high, (double)accuracy, &errorcode);

#else

          errorcode=0;
#ifdef QGGMRF
          UpdatedVoxelValue = CE_FunctionalSubstitution(low,high);

#else
          SurrogateUpdate=surrogateFunctionBasedMin();
          UpdatedVoxelValue=SurrogateUpdate;
#endif //QGGMRF
#endif//Surrogate function
          //printf("%lf\n",SurrogateUpdate);

          if(errorcode == 0)
          {
            //    printf("(%lf,%lf,%lf)\n",low,high,UpdatedVoxelValue);
            //  printf("Updated %lf\n",UpdatedVoxelValue);
#ifdef POSITIVITY_CONSTRAINT
            if(UpdatedVoxelValue < 0.0)
            { //Enforcing positivity constraints
              UpdatedVoxelValue = 0.0;
            }
#endif
          }
          else
          {
            if(THETA1 == 0 && low == 0 && high == 0) UpdatedVoxelValue = 0;
            else
            {
              printf("Error \n");
              printf("%d %d\n", j_new, k_new);
            }
          }

          //TODO Print appropriate error messages for other values of error code
          m_Geometry->Object->d[j_new][k_new][i] = UpdatedVoxelValue;
#ifdef NHICD
          MagUpdateMap->d[j_new][k_new] += fabs(m_Geometry->Object->d[j_new][k_new][i] - V);
#endif

#ifdef ROI
          if(Mask->d[j_new][k_new] == 1)
          {
            AverageUpdate += fabs(m_Geometry->Object->d[j_new][k_new][i] - V);
            AverageMagnitudeOfRecon += fabs(V); //computing the percentage update =(Change in mag/Initial magnitude)
          }
#endif

          //Update the ErrorSinogram
          for (uint32_t q = 0; q < TempMemBlock->count; q++)
          {
            uint16_t i_theta = floor(static_cast<float>(TempMemBlock->index[q] / (m_Sinogram->N_r)));
            uint16_t i_r = (TempMemBlock->index[q] % (m_Sinogram->N_r));
            VoxelLineAccessCounter = 0;
            for (uint32_t i_t = VoxelLineResponse[i].index[0]; i_t < VoxelLineResponse[i].index[0] + VoxelLineResponse[i].count; i_t++)
            //for(i_t = slice_index_min ; i_t <= slice_index_max; i_t++)
            {
              ErrorSino->d[i_theta][i_r][i_t] -= (NuisanceParams.I_0->d[i_theta]
                  * (TempMemBlock->values[q] * VoxelLineResponse[i].values[VoxelLineAccessCounter] * (m_Geometry->Object->d[j_new][k_new][i] - V)));
              VoxelLineAccessCounter++;
            }
          }
          Idx++;
        }
      }
      else
      {
        continue;
      }
    }
  }
  STOP_TIMER;
  PRINT_TIME;

#ifdef RANDOM_ORDER_UPDATES
  for (j = 0; j < m_Geometry->N_z; j++)
  { //Row index
    for (k = 0; k < m_Geometry->N_x; k++)
    { //Column index
      if(VisitCount->d[j][k] == 0)
      { printf("Pixel (%d %d) not visited\n", j, k);}}}
#endif

#ifdef COST_CALCULATE
  /*********************Cost Calculation***************************************************/
  cost.push_back(computeCost(ErrorSino, Weight));
  std::cout << indent << " cost[" << cost.size() << "] = " << cost.back() << std::endl;
  //printf("\n%lf", cost->d[cost_counter]);
  if(cost[cost.size()-1] - cost[cost.size() - 2] > 0)
  {
    std::cout << indent << "Cost just increased!" << std::endl;
    break;
  }
  fwrite( &(cost.back()), sizeof(DATA_TYPE), 1, Fp2);
  //  cost_counter++;
  /*******************************************************************************/
#else
  printf("%d\n",Iter);
#endif //Cost calculation endif
#ifdef ROI
  if(AverageMagnitudeOfRecon > 0)
  {
    printf("%d,%lf\n", Iter + 1, AverageUpdate / AverageMagnitudeOfRecon);
    if((AverageUpdate / AverageMagnitudeOfRecon) < m_Inputs->StopThreshold)
    {
      printf("This is the terminating point %d\n", Iter);
      m_Inputs->StopThreshold*=THRESHOLD_REDUCTION_FACTOR; //Reducing the thresold for subsequent iterations
      break;
    }
  }
#endif//ROI end
  if(getCancel() == true)
  {
    setErrorCondition(err);
    return;
  }

#ifdef WRITE_INTERMEDIATE_RESULTS

  if(Iter == NumOfWrites*WriteCount)
  {
    WriteCount++;
    sprintf(buffer,"%d",Iter);
    sprintf(Filename,"ReconstructedObjectAfterIter");
    strcat(Filename,buffer);
    strcat(Filename,".bin");
    Fp3 = fopen(Filename, "w");
    //  for (i=0; i < Geometry->N_y; i++)
    //    for (j=0; j < Geometry->N_z; j++)
    //      for (k=0; k < Geometry->N_x; k++)
    TempPointer = m_Geometry->Object;
    NumOfBytesWritten=fwrite(&(m_Geometry->Object->d[0][0][0]), sizeof(DATA_TYPE),m_Geometry->N_x*m_Geometry->N_y*m_Geometry->N_z, Fp3);
    printf("%d\n",NumOfBytesWritten);

    fclose(Fp3);
  }
#endif
  }
