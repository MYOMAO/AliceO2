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
#include <iostream>
#include <dirent.h> 
#include <stdio.h> 
#include <algorithm>
#include <iterator>
#include <chrono>
#include <thread>


using namespace o2::framework;
using namespace o2::itsmft;
using namespace std;

namespace o2
{
	namespace ITS
	{


		void RawPixelReader::init(InitContext& ic)
		{
			IndexPush = 0;

			std::ifstream EventPush("Config.dat");	
			EventPush >> EventPerPush >> TrackError >> workdir;


			LOG(INFO) << "EventPerPush = " << EventPerPush << "   TrackError = " << TrackError << "  workdir = " << workdir;


			o2::base::GeometryManager::loadGeometry();


			FolderNames = GetFName(workdir);

			cout << "NFolder = " << FolderNames.size() << endl;
			for (int i = 0; i < FolderNames.size(); i++){

				cout << "FDN = " << FolderNames[i] << endl;

				FileNames.push_back(GetFName(FolderNames[i]));

				cout << "FDN File Size = " << FileNames[i].size() << endl;


				for(int j = 0; j < FileNames[i].size(); j++){

					cout << "FDN File = " << FileNames[i][j] << endl;

				}

			}
			for(int i = 0; i < NError; i++){
				Error[i] = 0;
			}



			//		EventPerPush = 3000;
			EventRegistered = 0;
			TotalPixelSize = 0;

			//	GetFileName("infile");

			//
			o2::ITS::GeometryTGeo * geom = o2::ITS::GeometryTGeo::Instance ();
			geom->fillMatrixCache (o2::utils::bit2Mask (o2::TransformType::L2G));	
			const Int_t numOfChips = geom->getNumberOfChips ();	
			LOG(INFO) << "numOfChips = " << numOfChips;
			setNChips (numOfChips);	
			j = 0;
			FileDone = 1;
			PercentDone = 0.0;
			FileRemain = 0;







		}



		void RawPixelReader::run(ProcessingContext& pc)
		{


			cout << "------------------------------------------------------------------------------------------" << endl;
			cout << "------------------------------------------------------------------------------------------" << endl;

			cout << "New Cycle" << endl;	

			cout << "Wake Up Bro" << endl;

			cout << "Old Folder Size = " << FolderNames.size() << endl;


			NowFolderNames = GetFName(workdir);



			cout << "Now NFolder = " << NowFolderNames.size() << endl;
			for (int i = 0; i < NowFolderNames.size(); i++){

				//	cout << "FDN = " << NowFolderNames[i] << endl;

				NowFileNames.push_back(GetFName(NowFolderNames[i]));

				//cout << "Now FDN File Size = " << NowFileNames[i].size() << endl;


				for(int j = 0; j < NowFileNames[i].size(); j++){

					//	cout << "Now FDN File = " << NowFileNames[i][j] << endl;

				}

			}

			std::set_difference(NowFolderNames.begin(), NowFolderNames.end(), FolderNames.begin(), FolderNames.end(),std::inserter(DiffFolderName, DiffFolderName.begin()));

			cout << "Difference Size Between New and Initial Runs = " <<   DiffFolderName.size() << endl;


			if( DiffFolderName.size() == 0 ){
				cout << "No New Run -- No Need to Fucking Reset" << endl;
				ResetCommand = 0;
				pc.outputs().snapshot(Output{ "TST", "TEST", 0, Lifetime::Timeframe }, ResetCommand);
			}


			if( DiffFolderName.size() > 0){
				cout << "New Run Started -- Reset All Histograms" << endl;
				ResetCommand = 1;	
				pc.outputs().snapshot(Output{ "TST", "TEST", 0, Lifetime::Timeframe }, ResetCommand);
				for(int i = 0; i < NError; i++){
					Error[i] = 0;
				}
				ErrorVec.clear();
				ResetCommand = 0;	
				LOG(INFO) << "DONE Reset Histogram Decision";

			}


			LOG(INFO) << "Start Creating New Now Vector";



			LOG(INFO) << "Get IN LOOP";
			for(int i = 0;  i < FolderNames.size(); i++){
				std::set_difference(NowFileNames[i].begin(), NowFileNames[i].end(), FileNames[i].begin(), FileNames[i].end(),std::inserter(DiffFileNamePush, DiffFileNamePush.begin()));
				DiffFileNames.push_back(DiffFileNamePush);
				cout << "Difference File Size Between New and Initial Runs " <<   DiffFileNames[i].size() << endl;

				DiffFileNamePush.clear();
			}

			LOG(INFO) << "DONE GRABING Existing";

			for(int i = FolderNames.size();  i < NowFolderNames.size(); i++){
				DiffFileNames.push_back(NowFileNames[i]);
				cout << "New File Size Between New and Initial Runs " <<   DiffFileNames[i].size() << endl;
			}	

			LOG(INFO) << "Total New Files = " << DiffFileNames.size();

			LOG(INFO) << "DONE Creating Difference";			




			LOG(INFO) << "DiffFileNames Size = " << DiffFileNames.size();

			LOG(INFO) << "DONE Checking -- Reseting Vectors to the latest vector";



			LOG(INFO) << "DONE Updateing Vectors";


			LOG(INFO) << "Start Loop Bro";

			for (int i = 0; i < NowFolderNames.size(); i++){

				LOG(INFO) << "i = " << i << "    DiffFileNames[i].size() = " << DiffFileNames[i].size();			
				pos = NowFolderNames[i].find_last_of("/");

				if (pos != string::npos)   RunID =  NowFolderNames[i].substr(pos+1);


				//	for(int j = 0; j < DiffFileNames[i].size(); j++){
				if(DiffFileNames[i].size() > 0 && FileDone == 1){ 
					FileRemain = DiffFileNames[i].size();
					FileDone = 0;
					cout << "RunID = " << RunID << endl;
					cout << "File Location = " << DiffFileNames[i][0] << endl;


					size_t	last_index1 = RunID.find_last_not_of("0123456789");
					string RunIDS = RunID.substr(last_index1 + 1);


					string FileIDS;
					pos = DiffFileNames[i][0].find_last_of("/");
					if (pos != string::npos)   FileIDS = DiffFileNames[i][0].substr(pos+1);

					cout << "Before FileIDS = " << FileIDS << endl;

					size_t	last_index2 = FileIDS.find_last_not_of("0123456789");
					FileIDS =  FileIDS.substr(last_index2 + 1);

					RunName = std::stoi(RunIDS);
					FileID =  std::stoi(FileIDS);


					ofstream fout(Form("ErrorData/ErrorLogRun%d_File%d.dat",RunName,FileID));
					fout << " START OF ERROR REPORT For Run " <<  RunName << "  File " <<  FileID << endl;	

					inpName = DiffFileNames[i][0];


					//		pc.outputs().snapshot(Output{ "TST", "Run", 0, Lifetime::Timeframe }, RunID.data());
					//		pc.outputs().snapshot(Output{ "TST", "File", 0, Lifetime::Timeframe }, inpName.data());







					EventRegistered = 0;
					LOG(INFO) << "inpName = " << inpName;

					o2::itsmft::RawPixelReader<o2::itsmft::ChipMappingITS> rawReader;		
					rawReader.setPadding128(true); // payload GBT words are padded to 16B
					rawReader.setVerbosity(0);
					rawReader.setMinTriggersToCache(1025);

					rawReader.openInput(inpName);
					//mDigits.clear();
					//mMultiDigits.clear();

					int Index = 0;
					int IndexMax = -1;
					int NChip = 0;
					int NChipMax = -1;
					int TimePrint = 0;
					using RawReader=o2::itsmft::RawPixelReader<o2::itsmft::ChipMappingITS>;
					auto &rawErrorReader = reinterpret_cast<RawReader&>(rawReader);



					while (mChipData = rawReader.getNextChipData(mChips)) {
						if(NChip < NChipMax) break;
						//	cout << "Pass Chip" << endl;


						const auto* ruInfo = rawErrorReader.getCurrRUDecodeData()->ruInfo;	
						const auto& statRU =  rawErrorReader.getRUDecodingStatSW( ruInfo->idSW );

						const auto& pixels = mChipData->getData();
						int PixelSize = mChipData->getData().size();


						NEvent = statRU->nPackets;

						TotalPixelSize = TotalPixelSize + PixelSize;


						if(NEvent > (EventRegistered + 1) * EventPerPush){				
							NDigits.push_back(TotalPixelSize);
							EventRegistered = EventRegistered + 1;
							ErrorVec.push_back(Error);
							cout << "TotalPixelSize = " << TotalPixelSize << "  Pushed" << endl; 
							TotalPixelSize = 0;
						}

						if (NEvent%100000==0 && TimePrint == 0){
							cout << "Event Number = " << NEvent   << endl;
							TimePrint = 1;
						}

						if(NEvent%100000 != 0 ) TimePrint = 0;


						if(Error[0]  < 4294967295)  Error[0] = Error[0]  + (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrPageCounterDiscontinuity];
						if(Error[1]  < 4294967295)	Error[1] = Error[1]  + (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrRDHvsGBTHPageCnt];
						if(Error[2]  < 4294967295)	Error[2] = Error[2]  + (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrMissingGBTHeader];
						if(Error[3]  < 4294967295)	Error[3] = Error[3]  + (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrMissingGBTTrailer];
						if(Error[4]  < 4294967295)	Error[4] = Error[4]  + (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrNonZeroPageAfterStop];
						if(Error[5]  < 4294967295)	Error[5] = Error[5]  + (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrUnstoppedLanes];  
						if(Error[6]  < 4294967295)	Error[6] = Error[6]  + (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrDataForStoppedLane];
						if(Error[7]  < 4294967295)	Error[7] = Error[7]  + (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrNoDataForActiveLane];
						if(Error[8]  < 4294967295)	Error[8] = Error[8]  + (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrIBChipLaneMismatch];
						if(Error[9]  < 4294967295)	Error[9] = Error[9]  + (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrCableDataHeadWrong];
						if(Error[10]  < 4294967295)	Error[10] = Error[10]  + (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrPacketCounterJump];

						if(TrackError == 1){
							if(NEventPre != NEvent){
								if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrPageCounterDiscontinuity] > 0 || (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrRDHvsGBTHPageCnt] > 0 || (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrMissingGBTHeader] > 0 || (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrMissingGBTHeader] > 0||(int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrMissingGBTTrailer] > 0 || (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrNonZeroPageAfterStop] > 0 || (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrUnstoppedLanes] > 0 ||(int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrDataForStoppedLane] > 0 ||  (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrNoDataForActiveLane] > 0 ||  (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrIBChipLaneMismatch] > 0 || (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrIBChipLaneMismatch] > 0 || (int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrPacketCounterJump] > 0){
									fout << "Event Number = " << NEvent   << endl;
									fout << " ------------------------------------------------"   << endl;	
									if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrPageCounterDiscontinuity] > 0) fout << "Error ID 1: ErrPageCounterDiscontinuity Is Detected"  << endl;
									if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrRDHvsGBTHPageCnt] > 0) fout << "Error ID 2: ErrRDHvsGBTHPageCnt Is Detected"  << endl;
									if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrMissingGBTHeader] > 0) fout << "Error ID 3: ErrMissingGBTHeader Is Detected"  << endl;
									if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrMissingGBTTrailer] > 0) fout << "Error ID 4: ErrMissingGBTTrailer Is Detected"  << endl;
									if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrNonZeroPageAfterStop] > 0) fout << "Error ID 5: ErrNonZeroPageAfterStop Is Detected"  << endl;
									if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrUnstoppedLanes] > 0) fout << "Error ID 6: ErrUnstoppedLanes Is Detected"  << endl;
									if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrDataForStoppedLane] > 0) fout << "Error ID 7: ErrDataForStoppedLane Is Detected"  << endl;
									if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrNoDataForActiveLane] > 0) fout << "Error ID 8: ErrNoDataForActiveLane Is Detected"  << endl;
									if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrIBChipLaneMismatch] > 0) fout << "Error ID 9: ErrIBChipLaneMismatch Is Detected"  << endl;
									if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrIBChipLaneMismatch] > 0) fout << "Error ID 10: ErrCableDataHeadWrong Is Detected"  << endl;
									if((int)statRU->errorCounts[o2::itsmft::GBTLinkDecodingStat::ErrPacketCounterJump] > 0) fout << "Error ID 11: ErrPacketCounterJump Is Detected"  << endl;

									fout << " ------------------------------------------------"  << endl;
								}
							}
						}

						int ChipID = mChipData->getChipID();

						for (auto& pixel : pixels) {
							if(Index < IndexMax) break;
							//			cout << "Pass Pixel" << endl;	
							int col = pixel.getCol();
							int row = pixel.getRow();

							//			LOG(INFO) << "Chip ID Before " << ChipID << " Row = " << row << "   Column = " << col;

							mDigits.emplace_back(ChipID, NEvent, row, col, 0);
							//			LOG(INFO) << "Chip ID After " << mDigits[Index].getChipIndex() << " Row = " << mDigits[Index].getRow() << "   Column = " << mDigits[Index].getColumn();
							Index = Index + 1;
						}
						NChip = NChip + 1;
						NEventPre = NEvent;
					}
					cout << "Final TotalPixelSize = " << TotalPixelSize << endl; 
					NDigits.push_back(TotalPixelSize);
					ErrorVec.push_back(Error);
					LOG(INFO) <<"Run " << NowFolderNames[i] << " File " << inpName  <<  "    Integrated Raw Pixel Pushed " << mDigits.size();
					if(FolderNames.size() < NowFolderNames.size()) FileNames.push_back(NewNextFold);
					FileNames[i].push_back(inpName);
					fout << " END OF ERROR REPORT " << endl;	
				}
			}


			LOG(INFO) << "DONE Pushing";


			LOG(INFO) << "IndexPush Before = " << IndexPush << "  mDigits.size() =  "  << mDigits.size(); 
			
			if(mDigits.size() > 0) PercentDone = double(IndexPush)/double(mDigits.size());
			cout<< "Percentage Processed = " << Form("%.2f",100.*PercentDone) << endl;

	
			if(IndexPush < mDigits.size()){
				for(int i = 0; i < NDigits[j]; i++){
					mMultiDigits.push_back(mDigits[IndexPush + i]);
				}
				LOG(INFO) << "j = " << j << "   NDgits = " << NDigits[j]  << "    mMultiDigits Pushed = " << mMultiDigits.size();

				cout << "RunIDS = " << RunName << "   FileIDS = " << FileID << endl;

				pc.outputs().snapshot(Output{ "TST", "Run", 0, Lifetime::Timeframe }, RunName);
				pc.outputs().snapshot(Output{ "TST", "File", 0, Lifetime::Timeframe }, FileID);

				pc.outputs().snapshot(Output{ "TST", "Error", 0, Lifetime::Timeframe }, ErrorVec[j]);
				IndexPushEx =  IndexPush + NDigits[j];
				LOG(INFO) << "IndexPushEx = " << IndexPushEx << "  mDigits.size() " <<  mDigits.size();
				if(IndexPushEx >  mDigits.size() - 5) FileDone = 1;
				LOG(INFO) << "FileDone = " << FileDone;
				LOG(INFO) << "FileRemain = " << FileRemain;
	
				FileInfo = FileDone + FileRemain * 10;

				pc.outputs().snapshot(Output{ "TST", "Finish", 0, Lifetime::Timeframe }, FileInfo);

				/*
				   pc.outputs().snapshot(Output{ "TST", "Error0", 0, Lifetime::Timeframe }, ErrorVec[j][0]);
				   pc.outputs().snapshot(Output{ "TST", "Error1", 0, Lifetime::Timeframe }, ErrorVec[j][1]);
				   pc.outputs().snapshot(Output{ "TST", "Error2", 0, Lifetime::Timeframe }, ErrorVec[j][2]);
				   pc.outputs().snapshot(Output{ "TST", "Error3", 0, Lifetime::Timeframe }, ErrorVec[j][3]);
				   pc.outputs().snapshot(Output{ "TST", "Error4", 0, Lifetime::Timeframe }, ErrorVec[j][4]);
				   pc.outputs().snapshot(Output{ "TST", "Error5", 0, Lifetime::Timeframe }, ErrorVec[j][5]);
				   pc.outputs().snapshot(Output{ "TST", "Error6", 0, Lifetime::Timeframe }, ErrorVec[j][6]);
				   pc.outputs().snapshot(Output{ "TST", "Error7", 0, Lifetime::Timeframe }, ErrorVec[j][7]);
				   pc.outputs().snapshot(Output{ "TST", "Error8", 0, Lifetime::Timeframe }, ErrorVec[j][8]);
				   pc.outputs().snapshot(Output{ "TST", "Error9", 0, Lifetime::Timeframe }, ErrorVec[j][9]);
				   */

				//	LOG(INFO) << "mDigits.size() = " << mDigits.size();
				LOG(INFO) << "IndexPush = " << IndexPush << "    Chip ID Pushing " << mDigits[IndexPush].getChipIndex();
				//	pc.outputs().snapshot(Output{ "ITS", "DIGITS", 0, Lifetime::Timeframe }, mDigits[IndexPush++]);
				pc.outputs().snapshot(Output{ "ITS", "DIGITS", 0, Lifetime::Timeframe }, mMultiDigits);
				mMultiDigits.clear();
				IndexPush = IndexPush + NDigits[j];
				j = j + 1;
			}




			LOG(INFO) << "IndexPush After = " << IndexPush;






			/*
			   LOG(INFO) << "Before:  " << "IndexPush = " << IndexPush << "     mDigits.size() = " <<  mDigits.size(); 
			   while(IndexPush < mDigits.size()){
			//	LOG(INFO) << "mDigits.size() = " << mDigits.size();
			pc.outputs().snapshot(Output{ "ITS", "DIGITS", 0, Lifetime::Timeframe }, mDigits[IndexPush++]);
			if(IndexPush%100000==0) 	LOG(INFO) << "IndexPush = " << IndexPush << "    Chip ID Pushing " << mDigits[IndexPush].getChipIndex();
			}
			//pc.services().get<ControlService>().readyToQuit(true);
			LOG(INFO) << "After:  " << "IndexPush = " << IndexPush << "     mDigits.size() = " <<  mDigits.size(); 
			*/

			FolderNames.clear();
			//	FileNames.clear();
			NewNextFold.clear();
			FolderNames = NowFolderNames;
			//	FileNames = NowFileNames;

			NowFolderNames.clear();
			NowFileNames.clear();
			DiffFileNames.clear();
			DiffFolderName.clear();

			LOG(INFO) << "Pushing Reset Histogram Decision";

			cout << "Resetting Pushing Things" << endl;

			//Resetting a New File //
			if(IndexPush >  mDigits.size() - 5){
				mDigits.clear();
				IndexPush = 0;
				j = 0;
				NDigits.clear();
				FileDone = 1;
				pc.outputs().snapshot(Output{ "TST", "Finish", 0, Lifetime::Timeframe }, FileDone);	
				PercentDone = 0;
			}

			cout << "Start Sleeping Bro" << endl;
			cout << " " << endl;
			cout << " " << endl;
			cout << " " << endl;
			cout << " " << endl;
			cout << " " << endl;
			cout << " " << endl;
			cout << " " << endl;
			cout << " " << endl;
			cout << " " << endl;
			cout << " " << endl;	

			std::this_thread::sleep_for(std::chrono::milliseconds(1000));

			}



			std::vector<string> RawPixelReader::GetFName(std::string folder)
			{

				DIR           *dirp;
				struct dirent *directory;

				char cstr[folder.size()+1];
				strcpy(cstr, folder.c_str());
				dirp = opendir(cstr);
				std::vector<string> names;
				//string search_path = folder + "/*";
				if(dirp){

					while((directory = readdir(dirp)) != NULL){

						//printf("%s\n", directory->d_name);

						if ( !(!strcmp(directory->d_name, ".") || !strcmp(directory->d_name, ".."))) names.push_back(folder + "/" + directory->d_name);

					}

					closedir(dirp);
				}

				cout << "names size = " << names.size() << endl;
				return(names);
			}


			DataProcessorSpec getRawPixelReaderSpec()
			{
				return DataProcessorSpec{
					"Raw-Pixel-Reader",
						Inputs{},
						Outputs{
							OutputSpec{ "ITS", "DIGITS", 0, Lifetime::Timeframe },
							OutputSpec{ "TST", "TEST", 0, Lifetime::Timeframe },	
							OutputSpec{ "TST", "Error", 0, Lifetime::Timeframe },	
							OutputSpec{ "TST", "Run", 0, Lifetime::Timeframe },	
							OutputSpec{ "TST", "File", 0, Lifetime::Timeframe },
							OutputSpec{ "TST", "Finish", 0, Lifetime::Timeframe },		
							/*
							   OutputSpec{ "TST", "Error0", 0, Lifetime::Timeframe },
							   OutputSpec{ "TST", "Error1", 0, Lifetime::Timeframe },
							   OutputSpec{ "TST", "Error2", 0, Lifetime::Timeframe },
							   OutputSpec{ "TST", "Error3", 0, Lifetime::Timeframe },
							   OutputSpec{ "TST", "Error4", 0, Lifetime::Timeframe },
							   OutputSpec{ "TST", "Error5", 0, Lifetime::Timeframe },
							   OutputSpec{ "TST", "Error6", 0, Lifetime::Timeframe },
							   OutputSpec{ "TST", "Error7", 0, Lifetime::Timeframe },
							   OutputSpec{ "TST", "Error8", 0, Lifetime::Timeframe },
							   OutputSpec{ "TST", "Error9", 0, Lifetime::Timeframe },
							   */	
							//		OutputSpec{ "TST", "TEST3", 0, Lifetime::Timeframe },			
						},
						AlgorithmSpec{ adaptFromTask<RawPixelReader>() },
				};
			}



		}
	}
