#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <string>
#include "TTree.h"

#include "Framework/ControlService.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "ITSMFTReconstruction/GBTWord.h"
#include "ITSMFTReconstruction/PayLoadCont.h"
#include "ITSMFTReconstruction/PixelData.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITSMFTReconstruction/RawPixelReader.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "ITSRawWorkflow/RawPixelReaderSpec.h"
#include "DetectorsBase/GeometryManager.h"
#include <TCanvas.h>

using namespace o2::framework;
using namespace o2::ITSMFT;

namespace o2
{
	namespace ITS
	{


		void RawPixelReader::init(InitContext& ic)
		{
			IndexPush = 0;


			o2::base::GeometryManager::loadGeometry ();
			LOG(INFO) << "inpName = " << inpName;
			o2::ITS::GeometryTGeo * geom = o2::ITS::GeometryTGeo::Instance ();
			geom->fillMatrixCache (o2::utils::bit2Mask (o2::TransformType::L2G));	
			const Int_t numOfChips = geom->getNumberOfChips ();	
			LOG(INFO) << "numOfChips = " << numOfChips;
			setNChips (numOfChips);	
			rawReader.openInput(inpName);
			rawReader.setPadding128(true); // payload GBT words are padded to 16B
			//	rawReader.imposeMaxPage(1); // pages are 8kB in size (no skimming)
			rawReader.setVerbosity(0);
			mDigits.clear();
			mMultiDigits.clear();

			int Index = 0;
			int IndexMax = 10000;
			int NChip = 0;
			int NChipMax = 20;



			while (mChipData = rawReader.getNextChipData(mChips)) {
				if(NChip > NChipMax) break;

				const auto& pixels = mChipData->getData();
				int ChipID = mChipData->getChipID();

				for (auto& pixel : pixels) {
					if(Index > IndexMax) break;
					int col = pixel.getCol();
					int row = pixel.getRow();

					//			LOG(INFO) << "Chip ID Before " << ChipID << " Row = " << row << "   Column = " << col;

					mDigits.emplace_back(ChipID, Index, row, col, 0);
					//			LOG(INFO) << "Chip ID After " << mDigits[Index].getChipIndex() << " Row = " << mDigits[Index].getRow() << "   Column = " << mDigits[Index].getColumn();
					Index = Index + 1;


				}
				NChip = NChip + 1;

			}
			LOG(INFO) << "Integrated Raw Pixel Pushed " << mDigits.size();


		}



		void RawPixelReader::run(ProcessingContext& pc)
		{


			int NDigits = 10;

			if(IndexPush < mDigits.size()){

				for(int i = 0; i < NDigits; i++){

					mMultiDigits.push_back(mDigits[IndexPush + i]);
				}
				LOG(INFO) << "NDgits = " << NDigits  << "    mMultiDigits Pushed = " << mMultiDigits.size();

				//	LOG(INFO) << "mDigits.size() = " << mDigits.size();
				//	LOG(INFO) << "IndexPush = " << IndexPush << "    Chip ID Pushing " << mDigits[IndexPush].getChipIndex();
				//	pc.outputs().snapshot(Output{ "ITS", "DIGITS", 0, Lifetime::Timeframe }, mDigits[IndexPush++]);
				pc.outputs().snapshot(Output{ "ITS", "DIGITS", 0, Lifetime::Timeframe }, mMultiDigits);
				mMultiDigits.clear();
				IndexPush = IndexPush + NDigits;
			}





			/*

			   if(IndexPush < mDigits.size()){

			//	LOG(INFO) << "mDigits.size() = " << mDigits.size();
			//	LOG(INFO) << "IndexPush = " << IndexPush << "    Chip ID Pushing " << mDigits[IndexPush].getChipIndex();
			pc.outputs().snapshot(Output{ "ITS", "DIGITS", 0, Lifetime::Timeframe }, mDigits[IndexPush++]);
			}
			//pc.services().get<ControlService>().readyToQuit(true);
			*/

		}


		DataProcessorSpec getRawPixelReaderSpec()
		{
			return DataProcessorSpec{
				"Raw-Pixel-Reader",
					Inputs{},
					Outputs{
						OutputSpec{ "ITS", "DIGITS", 0, Lifetime::Timeframe },
					},
					AlgorithmSpec{ adaptFromTask<RawPixelReader>() },
			};
		}



	}
}
