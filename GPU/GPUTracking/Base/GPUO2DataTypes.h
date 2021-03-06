// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUO2DataTypes.h
/// \author David Rohr

#ifndef O2_GPU_GPUO2DATATYPES_H
#define O2_GPU_GPUO2DATATYPES_H

#ifdef HAVE_O2HEADERS
#include "DataFormatsTPC/ClusterNative.h"
#include "DetectorsBase/MatLayerCylSet.h"
#include "TRDBase/TRDGeometryFlat.h"
#else
namespace o2
{
namespace tpc
{
struct ClusterNative {
};
struct ClusterNativeAccess {
  const ClusterNative* clustersLinear;
  unsigned int nClusters[GPUCA_NSLICES][GPUCA_ROW_COUNT];
  unsigned int nClustersSector[GPUCA_NSLICES];
  unsigned int clusterOffset[GPUCA_NSLICES][GPUCA_ROW_COUNT];
  unsigned int nClustersTotal;
};
} // namespace tpc
namespace base
{
struct MatBudget {
};
class MatLayerCylSet
{
};
} // namespace base
namespace trd
{
class TRDGeometryFlat
{
};
} // namespace trd
} // namespace o2
#endif
#include "GPUdEdxInfo.h"

#endif
