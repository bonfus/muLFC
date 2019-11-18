#pragma once

struct CLattice; // An opaque type that we'll use as a handle
typedef struct CLattice CLattice;

typedef void(*CLatticeErrorHandlerFn)(const char * error_message, void * user_data);

typedef struct CLatticeError {
  CLatticeErrorHandlerFn eh;
  void * user_data;
} CLatticeError;

CLattice * lattice_create( const double *unitCell,
                           const int nAtoms,
                           const double *atomicPositions,
                           const double* in_K,
                           const double* in_fc,
                           const double* in_phi,
                           CLatticeError * eh  );

void lattice_set_occupations ( CLattice *h,
                           const double *atomicOccupations,
                           const int *atomicOccupationsGroups,
                           const double *sitesCorrelation,
                           CLatticeError * eh );

void lattice_destroy( CLattice * v, CLatticeError * eh );
