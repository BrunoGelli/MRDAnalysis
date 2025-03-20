
	#include <stdlib.h>     /* srand, rand */
	#include <time.h>       /* time */
	#include "TH1F.h"
	#include "TVirtualFFT.h"
	#include "TRandom.h"
	#include "TF1.h"
	#include "TTree.h"
	#include "TCanvas.h"
	#include "TStyle.h"
	#include "TMath.h"
	#include <math.h> 
	#include <stdio.h>
	#include <string.h>
	#include <TSpectrum.h>
	#include <TSpectrumTransform.h>

void graphtweaks();
void draw_trajectories(string name, int howMany, vector<double> x1, vector<double> y1, vector<double> z1, vector<double> x2, vector<double> y2, vector<double> z2);
double vector_elevation(double x, double y, double z);
double vector_azimuth(double x, double y, double z);
double vector_with_normal(double x, double y, double z);
void Normalize2DHistogram(TH2D* h2);
void draw_simple_cut_histograms(string name, string variable, vector<double> data_pre_cut, vector<double> data_post_cut, bool invert, int bin, int min, int max, TPaveLabel* label);
void cutMinMax(string CutName, vector<bool>& cut, vector<double> measurement, double min, double max, TPaveLabel* label);
void cutString(string CutName, vector<bool>& cut, vector<string> measurement, string condition, TPaveLabel* label);
void cutValue(string CutName,  vector<bool>& cut, vector<int> measurement, int value, TPaveLabel* label);

struct ChannelInfo {
	int DetSys, Orient, Layer, Side, Number, Rack, Slot, Ch;
};

ChannelInfo getChannelInfo(const map<int, ChannelInfo>& channelMap, int channel);
map<int, ChannelInfo> readChannelMapping(const string& filename);
vector<TH2D*> convertHistogram(TH2D* originalHist, const map<int, ChannelInfo>& channelMap, int rack);
void drawFancyEffMap(TH2D* originalHist, int rack);
// void drawFitParameters(TF1 * fitFunction);

void MrdAnalysis()
{

	graphtweaks();

 	std::vector<TString> fileNames;
 	fileNames.push_back("mrd_observed_new_R3347S0_Beam.root");

    std::vector<double> Elevation,  Azimuth,  Normal;
    std::vector<double> sElevation, sAzimuth, sNormal;

 	std::vector<TVector3> directionVector;

	std::vector<double> 			StartVertexX, StartVertexY, StartVertexZ;	//declaring the hit position vectors
	std::vector<double> 			StopVertexX,  StopVertexY,  StopVertexZ;	//declaring the hit position vectors
	std::vector<double> 			sXs, sYs, sZs, sXf, sYf, sZf;				//declaring the hit position vectors - post cut

	std::vector<int> 				RunNum, 			sRunNum;
	std::vector<int> 				EvNum, 				sEvNum;
	std::vector<string> 			MRDTriggerType, 	sMRDTriggerType;
	std::vector<int> 				NumMrdTracks, 		sNumMrdTracks;
	std::vector<int> 				MrdTrackID, 		sMrdTrackID;
	std::vector<int> 				NPMTsHit, 			sNPMTsHit;
	std::vector<int> 				NLayersHit, 		sNLayersHit;
	std::vector<int> 				IsThrough, 			sIsThrough;
	std::vector<std::vector<int>> 	MissingChannel, 	sMissingChannel;
	std::vector<std::vector<int>> 	ExpectedChannel, 	sExpectedChannel;
	std::vector<int> 				NumMissingChannel, 	sNumMissingChannel;
	std::vector<std::vector<int>> 	MissingLayer, 		sMissingLayer;
	std::vector<int> 				NumMissingLayer, 	sNumMissingLayer;

	std::vector<bool> PassedCut;
	
	double			auxStartVertexX, auxStartVertexY, auxStartVertexZ;
	double			auxStopVertexX,  auxStopVertexY,  auxStopVertexZ;

	int 			auxRunNum;
	int				auxEvNum;
	string		   *auxMRDTriggerType	= nullptr;
	int				auxNumMrdTracks;
	int				auxMrdTrackID;
	int				auxNPMTsHit;
	int				auxNLayersHit;
	int				auxIsThrough;
	vector<int>	   *auxMissingChannel 	= nullptr;
	vector<int>	   *auxExpectedChannel 	= nullptr;
	int				auxNumMissingChannel;
	vector<int>	   *auxMissingLayer 	= nullptr;
	int				auxNumMissingLayer;

	TFile *f;
	TTree *T;

 	for (int i = 0; i < fileNames.size(); ++i)
 	{
 		
		f = new TFile(fileNames[i],"READ");
		
		if ( f->IsOpen() )
		{
			printf("File opened successfully\n");
		} 
		else 
		{
			printf("Error! Invalid File Name\n");
			return;
		}

		T = (TTree*)f -> Get("tree_trackfit");

		T->SetBranchAddress(            "RunNum",            &auxRunNum);
		T->SetBranchAddress(             "EvNum",             &auxEvNum);
		T->SetBranchAddress(    "MRDTriggerType",    &auxMRDTriggerType);
		T->SetBranchAddress(      "NumMrdTracks",      &auxNumMrdTracks);
		T->SetBranchAddress(        "MrdTrackID",        &auxMrdTrackID);
		T->SetBranchAddress(          "NPMTsHit",          &auxNPMTsHit);
		T->SetBranchAddress(        "NLayersHit",        &auxNLayersHit);
		T->SetBranchAddress(      "StartVertexX",      &auxStartVertexX);
		T->SetBranchAddress(      "StartVertexY",      &auxStartVertexY);
		T->SetBranchAddress(      "StartVertexZ",      &auxStartVertexZ);
		T->SetBranchAddress(       "StopVertexX",       &auxStopVertexX);
		T->SetBranchAddress(       "StopVertexY",       &auxStopVertexY);
		T->SetBranchAddress(       "StopVertexZ",       &auxStopVertexZ);
		T->SetBranchAddress(         "IsThrough",         &auxIsThrough);
		T->SetBranchAddress(    "MissingChannel",    &auxMissingChannel);
		T->SetBranchAddress(   "ExpectedChannel",   &auxExpectedChannel);
		T->SetBranchAddress( "NumMissingChannel", &auxNumMissingChannel);
		T->SetBranchAddress(      "MissingLayer",      &auxMissingLayer);
		T->SetBranchAddress(   "NumMissingLayer",   &auxNumMissingLayer);
				
		int nentries = (int)T->GetEntries();

		for (int i = 0; i < nentries; ++i)
		{
			T->GetEntry(i);

			RunNum.push_back(auxRunNum);
			EvNum.push_back(auxEvNum);
			MRDTriggerType.push_back(*auxMRDTriggerType);
			NumMrdTracks.push_back(auxNumMrdTracks);
			MrdTrackID.push_back(auxMrdTrackID);	
			NPMTsHit.push_back(auxNPMTsHit);
			NLayersHit.push_back(auxNLayersHit);
			StartVertexX.push_back(auxStartVertexX);
			StartVertexY.push_back(auxStartVertexY);
			StartVertexZ.push_back(auxStartVertexZ);
			StopVertexX.push_back(auxStopVertexX);
			StopVertexY.push_back(auxStopVertexY);
			StopVertexZ.push_back(auxStopVertexZ);
			IsThrough.push_back(auxIsThrough);
			MissingChannel.push_back(*auxMissingChannel);
			ExpectedChannel.push_back(*auxExpectedChannel);
			NumMissingChannel.push_back(auxNumMissingChannel);
			MissingLayer.push_back(*auxMissingLayer);
			NumMissingLayer.push_back(auxNumMissingLayer);

		}

 	}


 	cout << "There were " << EvNum.size() << " events in this file." << endl;

 	for (int i = 0; i < EvNum.size(); ++i)
 	{
 		PassedCut.push_back(true);
 	}

 	for (int i = 0; i < EvNum.size(); ++i)
 	{
		TVector3 start(StartVertexX[i], StartVertexY[i], StartVertexZ[i]); 
    	TVector3 end(   StopVertexX[i],  StopVertexY[i],  StopVertexZ[i]); 
    	directionVector.push_back(end - start);
 	}

    for (int i = 0; i < EvNum.size(); ++i)
    {
		Elevation.push_back(vector_elevation(directionVector[i].X(),directionVector[i].Y(),directionVector[i].Z()));
		Azimuth.push_back(vector_azimuth(directionVector[i].X(),directionVector[i].Y(),directionVector[i].Z()));
		Normal.push_back(vector_with_normal(directionVector[i].X(),directionVector[i].Y(),directionVector[i].Z()));
    }

    TPaveLabel* paveLabel = new TPaveLabel(0.15,0.70,0.35,0.90, "Cut Conditions:", "NDC");
    paveLabel->SetFillColor(0);   // Transparent background
    // paveLabel->SetTextAlign(12);  // Align text to top-left
    paveLabel->SetBorderSize(1);  // Thin border
	paveLabel->SetTextSize(0.1); // text size

    cutString("Event Type Selector",   PassedCut, MRDTriggerType, "Cosmic", paveLabel);
    cutValue("Throughgoing Selector",  PassedCut, IsThrough,             1, paveLabel);
    cutMinMax("Azimuth Selector",  	   PassedCut, Azimuth,         -10, 10, paveLabel);
    // cutMinMax("Elevation Selector",  	   PassedCut, Elevation,   -15, 15, paveLabel);

    for (int i = 0; i < EvNum.size(); ++i)
    {
    	if (PassedCut[i])
    	{
    		sRunNum.push_back(RunNum[i]); 
    		sEvNum.push_back(EvNum[i]); 
    		sMRDTriggerType.push_back(MRDTriggerType[i]);
			sNumMrdTracks.push_back(NumMrdTracks[i]); 
			sMrdTrackID.push_back(MrdTrackID[i]);
			sNPMTsHit.push_back(NPMTsHit[i]); 
			sNLayersHit.push_back(NLayersHit[i]); 
			sIsThrough.push_back(IsThrough[i]);
			sMissingChannel.push_back(MissingChannel[i]); 
			sExpectedChannel.push_back(ExpectedChannel[i]);	
			sNumMissingChannel.push_back(NumMissingChannel[i]);
			sMissingLayer.push_back(MissingLayer[i]);     
			sNumMissingLayer.push_back(NumMissingLayer[i]);
    		sXs.push_back(StartVertexX[i]);	
    		sYs.push_back(StartVertexY[i]);	
    		sZs.push_back(StartVertexZ[i]);
    		sXf.push_back(StopVertexX[i]);	
    		sYf.push_back(StopVertexY[i]);	
    		sZf.push_back(StopVertexZ[i]);
    		sElevation.push_back(Elevation[i]); 
    		sAzimuth.push_back(Azimuth[i]); 
    		sNormal.push_back(Normal[i]);
    	}
    }

//************ Output File ************	

	TFile outFile("outputFile.root", "RECREATE");

//************ First ploting the number a number N of events in X,Y,Z ************	

	string drawPreCut = "preCut", drawPostCut = "postCut";

	draw_trajectories(drawPreCut,  1000, StartVertexX, StartVertexY, StartVertexZ, StopVertexX, StopVertexY, StopVertexZ);

	draw_trajectories(drawPostCut, 1000, sXs, sYs, sZs, sXf, sYf, sZf);

//************ Now lets plot the major variables for this analysis - before and after the cut  ************

	TCanvas *Simple = new TCanvas("Simple", "Simple plots", 4096, 2160);
	Simple->Divide(2,2);

	Simple->cd(1);
	draw_simple_cut_histograms("Elevation Forward Going", "Elevation", Elevation, sElevation, false, 45, 0, 90, paveLabel);
	
	Simple->cd(2);
	draw_simple_cut_histograms("Elevation Backward Going", "Elevation", Elevation, sElevation, true, 45, 0, 90, paveLabel);
	
	Simple->cd(3);
	draw_simple_cut_histograms("Azimuth", "Azimuth", Azimuth, sAzimuth, false, 90, -90, 90, paveLabel);
	
	Simple->cd(4);
	draw_simple_cut_histograms("Normal",  "Normal", Normal, sNormal, false, 45, 0, 90, paveLabel);

	Simple->Write("Simple Geometry Histograms");
	Simple->SaveAs("SimpleGeometryHistograms.png");

//************ Now lets plot a 2D Histogram of the missed channels  ************

	TH2D* missedChannel = new TH2D("Missed Chan", "Missed Chan", 18, 0, 90, 270, 40, 310);

	for (Int_t i=0; i<sMissingChannel.size(); i++)
	{

		if (sNormal[i] > 0)
		{
			for (int j = 0; j < sMissingChannel[i].size(); ++j)
			{
				missedChannel->Fill(sNormal[i], sMissingChannel[i][j]);
			}
		}
		if (sNormal[i] < 0)
		{
			for (int j = 0; j < sMissingChannel[i].size(); ++j)
			{
				missedChannel->Fill(-sNormal[i], sMissingChannel[i][j]);
			}
		}
	}

	missedChannel->SetTitle("Missed Channel Hits");
	missedChannel->GetXaxis()->SetTitle("Angle with Normal (#circ)");
	missedChannel->GetYaxis()->SetTitle("Channel Number");
	missedChannel->SetLineWidth(2);
	// missedChannel->SetStats(0);


	TCanvas *missedChannelCanvas = new TCanvas("missedChannel", "missedChannel", 4096, 2160);

	missedChannelCanvas->cd();
	missedChannel->Draw("lego2");
	paveLabel->Draw();

//************ Now lets plot a 2D Histogram of the expected channels  ************

	TH2D* expectedChannel = new TH2D("expected Chan", "expected Chan", 18, 0, 90, 270, 40, 310);

	for (Int_t i=0; i<sExpectedChannel.size(); i++)
	{

		if (sNormal[i] > 0)
		{
			for (int j = 0; j < sExpectedChannel[i].size(); ++j)
			{
				expectedChannel->Fill(sNormal[i], sExpectedChannel[i][j]);
			}
		}
		if (sNormal[i] < 0)
		{
			for (int j = 0; j < sExpectedChannel[i].size(); ++j)
			{
				expectedChannel->Fill(-sNormal[i], sExpectedChannel[i][j]);
			}
		}
	}

	expectedChannel->SetTitle("Expected Channel Hits");
	expectedChannel->GetXaxis()->SetTitle("Angle with Normal (#circ)");
	expectedChannel->GetYaxis()->SetTitle("Channel Number");
	expectedChannel->SetLineWidth(2);
	// expectedChannel->SetStats(0);


	TCanvas *expectedChannelCanvas = new TCanvas("expectedChannel", "expectedChannel", 4096, 2160);

	expectedChannelCanvas->cd();
	expectedChannel->Draw("lego2");
	paveLabel->Draw();

//************ Finally, using the 2 previous 2D Histograms to find the efficiency  ************

	// Create a new histogram for efficiency
	TH2D* hEfficiency = (TH2D*)missedChannel->Clone("Efficiency");
	hEfficiency->SetTitle("Efficiency Map; Normal angle (#circ); Channel Number; Efficiency");

	// Loop over bins to compute efficiency
	for (int i = 1; i <= hEfficiency->GetNbinsX(); ++i) {
		for (int j = 1; j <= hEfficiency->GetNbinsY(); ++j) {
		    double total = expectedChannel->GetBinContent(i, j);
		    double misses = missedChannel->GetBinContent(i, j);

		    // Efficiency calculation with protection against division by zero 
		    // also making sure we have some stat
		    if (total > 10)  //<---------------- change here to move the min stats threshold
		    { 
		        hEfficiency->SetBinContent(i, j, (total - misses) / total);
		    } 
		    else
		    {
		        hEfficiency->SetBinContent(i, j, 0);  // Or NAN if you prefer to mark invalid entries
		    }
		}
    }

	hEfficiency->SetStats(0);

    TCanvas* EffCanvas = new TCanvas("EfficiencyHist", "EfficiencyHist", 4096, 2160);
    
	gPad->SetRightMargin(.15);	

    hEfficiency->Draw("COLZ");
	// paveLabel->DrawPaveLabel(0.55,0.70,0.75,0.90, "Cut Conditions:", "");

	hEfficiency->Write();
	EffCanvas->Write();
	EffCanvas->SaveAs("EfficiencyHist.png");

//************ Lastly, lets summirize this into a #number of paddes missed vs angle 2D Histogram ************


	TCanvas *C_PaddleAnalysis2 = new TCanvas("C_PaddleAnalysis2", "Parede frontal",0,0,1920,1080);
	C_PaddleAnalysis2->cd();



	gPad->SetRightMargin(.15);	
	TH2D* PaddleAnalysis2 = new TH2D("", "", 12, 0, 60, 16., 0, 0.8);

	for (int i = 0; i < sNormal.size(); ++i)
	{
		if (sExpectedChannel[i].size()!=0)
		{
			double missed = sMissingChannel[i].size();
			double total  = sExpectedChannel[i].size();
			double percentMissed = (missed)/total;
			PaddleAnalysis2->Fill(sNormal[i], percentMissed);
		}
	}

    auto graph = new TGraphErrors(PaddleAnalysis2->GetNbinsX());

    vector<double> xVal, yVal, yUnc;

    // Compute mean and standard deviation for each x-bin
    for (int i = 1; i <= PaddleAnalysis2->GetNbinsX(); ++i) 
    {
        double xCenter = PaddleAnalysis2->GetXaxis()->GetBinCenter(i);

        // Collect all Y values in this X-bin
        std::vector<double> y_values;

        for (int j = 1; j <= PaddleAnalysis2->GetNbinsY(); ++j) 
        {
        	int bin = PaddleAnalysis2->GetBin(i, j);
            double entries = PaddleAnalysis2->GetBinContent(bin);


            for (int k = 0; k < entries; ++k) 
            {
                double y = PaddleAnalysis2->GetYaxis()->GetBinCenter(j);
                y_values.push_back(y);
            }
        }

        // Calculate mean and standard deviation
        double mean = 0.0, stddev = 0.0;
        if (!y_values.empty()) 
        {
            mean = TMath::Mean(y_values.begin(), y_values.end());
            stddev = TMath::StdDev(y_values.begin(), y_values.end());

        }

        graph->SetPoint(i - 1, xCenter, mean);
        graph->SetPointError(i - 1, 0, stddev);
    }

    // Fit with a linear function
    TF1 *fitFunc = new TF1("fitFunc", "[0]*e^(x*[1])+[2]", 0, 63);
    graph->Fit(fitFunc, "R");  // "R" restricts fit to graph range
    Normalize2DHistogram(PaddleAnalysis2);

    // Drawing
    PaddleAnalysis2->GetXaxis()->SetTitle("Angle with Normal (#circ)");
	PaddleAnalysis2->GetYaxis()->SetTitle("Fraction of missed paddles");
	PaddleAnalysis2->GetZaxis()->SetTitle("Normalized frequency (by the angle)");
	PaddleAnalysis2->SetStats(kFALSE);
    PaddleAnalysis2->Draw("COLZ");          // Draw 2D histogram as background
    graph->SetMarkerStyle(20); 
    graph->SetMarkerColor(1);
    graph->Draw("P SAME");     // Draw points on top
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("SAME");     // Fitted line
    // drawFitParameters(fitFunc);

	PaddleAnalysis2->Write();
    C_PaddleAnalysis2->Write();
	C_PaddleAnalysis2->SaveAs("IntegratedPaddleResult.png");

//************ Final thing is to plot it vs the channel mapping ************

	map<int, ChannelInfo> channelMap = readChannelMapping("mapping.csv");


	vector<TH2D*> AngleDependentHistRack7 = convertHistogram(hEfficiency, channelMap, 7);
	vector<TH2D*> AngleDependentHistRack8 = convertHistogram(hEfficiency, channelMap, 8);

		for (auto hist : AngleDependentHistRack7) 
		{
		    hist->Write();
		    drawFancyEffMap(hist, 7);
		}
		for (auto hist : AngleDependentHistRack8) 
		{
		    hist->Write();
		    drawFancyEffMap(hist, 8);
		}


		outFile.Close();
		return;

} 

//************ ************************************************ ************
//************ 					Functions 						************
//************ ************************************************ ************

void graphtweaks()
{
	gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.95);
	gStyle->SetStatW(0.3);
	gStyle->SetStatH(0.2);
	gStyle->SetStatBorderSize(3);
	gStyle->SetLineColor(1);
	// gStyle->SetPalette(62); 
}

void draw_trajectories(string name, int howMany, vector<double> x1, vector<double> y1, vector<double> z1, vector<double> x2, vector<double> y2, vector<double> z2) 
{


    // First we get N random events between 0 and x1.size()
 	int vectorSize = x1.size();
 	vector<int> selectedEvents;

	TRandom3 *randomGenerator = new TRandom3();
	randomGenerator->SetSeed(time(0));

 	for (int i = 0; i < howMany; ++i)
 	{
 		selectedEvents.push_back(randomGenerator->Integer(vectorSize));
 	}
    int first_event = 0;
    int  last_event = 10;


	//now for plotting

    TCanvas *c1 = new TCanvas(name.c_str(),  name.c_str(), 800, 600);
    
    // Create an empty 3D graph (for better visualization)
    TGraph2D *Points = new TGraph2D();
    TGraph2D *arrowHead = new TGraph2D();

    Points->SetTitle("Track Vectors;X;Y;Z");

    for (int j = 0; j < selectedEvents.size(); ++j)
    {
    	int i = selectedEvents[j];
        Points->SetPoint(Points->GetN(),  x1[i], y1[i], z1[i]);
        Points->SetPoint(Points->GetN(),  x2[i], y2[i], z2[i]);
        arrowHead->SetPoint(arrowHead->GetN(),   x2[i], y2[i], z2[i]);

    }


    // Draw the 3D graph and lines
    Points->SetMarkerStyle(20);
    Points->SetMarkerSize(1);
    Points->Draw("P"); // "P" draws points

	arrowHead->SetMarkerStyle(29);  // Filled star (looks like an arrow tip)
	arrowHead->SetMarkerSize(2.0);
	arrowHead->SetMarkerColor(kRed);
	arrowHead->Draw("P same");

	// now let draw the lines
    for (int j = 0; j < selectedEvents.size(); ++j)
    {
    	int i = selectedEvents[j];
    	// Use TPolyLine3D to draw a line from start to end
        TPolyLine3D *line = new TPolyLine3D(2);
        line->SetPoint(0, x1[i], y1[i], z1[i]); // Start point
        line->SetPoint(1, x2[i], y2[i], z2[i]);       // End point
        line->SetLineColor(kRed);
        line->SetLineWidth(1);
        line->Draw("a1 same");
    }

    c1->Update();
    c1->Write(name.c_str());
}

void cutMinMax(string CutName, vector<bool>& cut, vector<double> measurement, double min, double max, TPaveLabel* label)
{
	if (cut.size() != measurement.size())
	{
		cout << "something is going wrong in the data (worst error text possible)" << endl;
	}

	int NumCut = 0;
	for (int i = 0; i < measurement.size(); ++i)
	{
		if (cut[i])
		{
			if (measurement[i] < min || measurement[i] > max)
			{
				cut[i] = false;
				NumCut++;
			}
		}
	}

	cout << "Number of events removed by " << CutName << ": "<< NumCut  << endl;
	string spacer = "         ";
	string legend = CutName + ": ";
	string legendLine2 = spacer + "Min = " + Form("%.0f", min) + ";";
	string legendLine3 = spacer + "Max = " + Form("%.0f", max) + ";";

    // Append the cut condition using #splitline{}
    if (label) {
        string currentText = label->GetLabel();

        if (currentText.empty() || currentText.find("#splitline") == string::npos) {
            // Start fresh or initialize first entry
            currentText = "#splitline{" + CutName + ":  Min = " + min + " | Max = "+ max + "}{}";
        } else {
            // Add new cuts with nested #splitline
            currentText = "#splitline{" + currentText + "}{" + CutName + ":  Min = " + min + " | Max = "+ max + "}";
        }

        label->SetLabel(currentText.c_str());
    }
}

void cutString(string CutName, vector<bool>& cut, vector<string> measurement, string condition, TPaveLabel* label)
{
	if (cut.size() != measurement.size())
	{
		cout << "something is going wrong in the data (worst error text possible)" << endl;
	}

	int NumCut = 0;

	for (int i = 0; i < measurement.size(); ++i)
	{
		if (cut[i])
		{
			if (measurement[i] != condition)
			{
				cut[i] = false;
				NumCut++;
			}
		}
	}

	cout << "Number of events removed by " << CutName << ": "<< NumCut  << endl;

    // Append the cut condition using #splitline{}
    if (label) {
        string currentText = label->GetLabel();

        if (currentText.empty() || currentText.find("#splitline") == string::npos) {
            // Start fresh or initialize first entry
            currentText = "#splitline{" + CutName + ": " + condition + "}{}";
        } else {
            // Add new cuts with nested #splitline
            currentText = "#splitline{" + currentText + "}{" + CutName + ": " + condition + "}";
        }

        label->SetLabel(currentText.c_str());
    }
}

void cutValue(string CutName, vector<bool>& cut, vector<int> measurement, int value, TPaveLabel* label)
{
	if (cut.size() != measurement.size())
	{
		cout << "something is going wrong in the data (worst error text possible)" << endl;
	}

	int NumCut = 0;
	for (int i = 0; i < measurement.size(); ++i)
	{
		if (cut[i])
		{
			if (measurement[i] != value)
			{
				cut[i] = false;
				NumCut++;
			}
		}
	}

	cout << "Number of events removed by " << CutName << ": "<< NumCut  << endl;


	string legend = CutName + ": " + Form("%d", value) + ";";

   if (label) {
        string currentText = label->GetLabel();

        if (currentText.empty() || currentText.find("#splitline") == string::npos) {
            // Start fresh or initialize first entry
            currentText = "#splitline{" + CutName + ": " + value + "}{}";
        } else {
            // Add new cuts with nested #splitline
            currentText = "#splitline{" + currentText + "}{" + CutName + ": " + value + "}";
        }

        label->SetLabel(currentText.c_str());
    }
}

double vector_elevation(double x, double y, double z) 
{
    double theta = atan(y/z) * 180.0 / 3.1415;  // Inclination angle
	return theta;
}

double vector_azimuth(double x, double y, double z) 
{
    double phi = atan(x/z) * 180.0 / 3.1415;  // Azimuthal angle in XZ-plane
	return phi;
}

double vector_with_normal(double x, double y, double z) 
{
	TVector3 Normal(0, 0, 1);  	// normal to the MRD
    TVector3 Trajectory(x, y, z);  	// Trajectory vector

 	// Calculate the dot product
    double dotProduct = Normal.Dot(Trajectory);

    // Calculate the magnitudes of A and B
    double magNormal = Normal.Mag();
    double magTrajectory = Trajectory.Mag();

    // Calculate the cosine of the angle
    double cosTheta = dotProduct / (magNormal * magTrajectory);

    // Ensure the value is within the range [-1, 1] due to floating-point precision
    cosTheta = std::max(-1.0, std::min(1.0, cosTheta));

    // Calculate the angle
    double angle = acos(cosTheta)* 180.0 / 3.1415;

	return angle;
}


void draw_simple_cut_histograms(string name, string variable, vector<double> data_pre_cut, vector<double> data_post_cut, bool invert, int bin, int min, int max, TPaveLabel* label)
{

	string cutName = name + " after cuts";
	string xAxName = variable +  "(#circ)";
	string cutVar  = variable +  " after cuts";

	TH1D *preCutHist = new TH1D(name.c_str(),	name.c_str(), 		bin,   min,  max);
	TH1D *posCutHist = new TH1D(cutName.c_str(),	cutName.c_str(), 	bin,   min,  max);
	

	for (Int_t i=0; i<data_pre_cut.size(); i++)
	{
		if (invert)
		{
			preCutHist->Fill(-data_pre_cut[i]);
		}
		else
		{
			preCutHist->Fill(data_pre_cut[i]);
		}
		
	}

	for (Int_t i=0; i<data_post_cut.size(); i++)
	{
		if (invert)
		{
			posCutHist->Fill(-data_post_cut[i]);
		}
		else
		{
			posCutHist->Fill(data_post_cut[i]);
		}
		
	}

	preCutHist->SetTitle(name.c_str());
	preCutHist->GetXaxis()->SetTitle(variable.c_str());
	preCutHist->GetYaxis()->SetTitle("#");
	preCutHist->SetLineWidth(2);
	preCutHist->Draw();
	preCutHist->SetStats(0);
	posCutHist->SetLineWidth(2);
	posCutHist->SetLineColor(2);
	posCutHist->Draw("same");

	auto Legend = new TLegend(0.65,0.65,0.85,0.85);
	Legend->AddEntry(preCutHist, variable.c_str() , "l");
	Legend->AddEntry(preCutHist, Form("Entries: %.0f", preCutHist->GetEntries()), "");

	Legend->AddEntry(posCutHist, cutVar.c_str() , "l");
	Legend->AddEntry(posCutHist, Form("Entries: %.0f", posCutHist->GetEntries()), "");

	Legend->Draw();
	label->Draw();

	preCutHist->Write(name.c_str());
	posCutHist->Write(cutName.c_str());

}

void Normalize2DHistogram(TH2D* h2) {
    if (!h2) {
        std::cerr << "Error: Histogram pointer is null!" << std::endl;
        return;
    }

    int nBinsX = h2->GetNbinsX(); // Number of bins along the x-axis
    int nBinsY = h2->GetNbinsY(); // Number of bins along the y-axis

    // Loop over x-axis bins (angle of incidence)
    for (int i = 1; i <= nBinsX; ++i) {
        double columnSum = 0.0;

        // First pass: Calculate total entries in the column
        for (int j = 1; j <= nBinsY; ++j) {
            columnSum += h2->GetBinContent(i, j);
        }

        // Second pass: Normalize each y-bin value
        if (columnSum > 0) {
            for (int j = 1; j <= nBinsY; ++j) {
                double currentValue = h2->GetBinContent(i, j);
                h2->SetBinContent(i, j, currentValue / columnSum);
            }
        }
    }
}



// Function to read the CSV and store channel mappings
map<int, ChannelInfo> readChannelMapping(const string& filename) 
{
    map<int, ChannelInfo> channelMap;
    ifstream file(filename);

    if (!file.is_open()) 
    {
        cerr << "Error: Could not open file " << filename << endl;
        return channelMap;
    }

    string line;
    // Skip header
    // getline(file, line);
    // getline(file, line);
    // getline(file, line);
    // getline(file, line);
    // getline(file, line);
    // getline(file, line);
    // getline(file, line);
    // getline(file, line);

    // Read data
    while (getline(file, line)) 
    {
        stringstream ss(line);
        vector<string> row;
        string value;

        // Split CSV line into individual values
        while (getline(ss, value, ',')) 
        {
            row.push_back(value);
        }

        if (row.size() >= 16)
        {  // Ensure row has enough columns
            int ChKey 	= stoi(row[0]);
            int DetSys 	= stoi(row[2]);		// 0 = FMV and 1 = MRD
            int Orient 	= stoi(row[3]);		// 0 = Horizontal and 1 = Vertical
            int Layer 	= stoi(row[4]);
            int Side 	= stoi(row[5]);		// 0=left or 1=right  || 0=bottom or 1=top
            int Number 	= stoi(row[6]);		// from left to right || from bottom to top	
            int Rack 	= stoi(row[13]);  	// Rack 
            int Slot 	= stoi(row[14]);  	// Slot Column 
            int Ch 		= stoi(row[15]);    // Channel Column
            channelMap[ChKey] = {DetSys, Orient, Layer, Side, Number, Rack, Slot, Ch};
        }
    }

    file.close();
    return channelMap;
}

// Function to retrieve Slot and Ch for a given channel
ChannelInfo getChannelInfo(const map<int, ChannelInfo>& channelMap, int channel) 
{
    if (channelMap.find(channel) != channelMap.end()) 
    {
        return channelMap.at(channel);
    } 
    else 
    {
        cerr << "Error: Channel " << channel << " not found in mapping." << endl;
        return {-1, -1, -1, -1, -1, -1, -1, -1}; // Return invalid values
    }
}

// Function to convert TH2D histogram into multiple Slot vs Ch histograms
vector<TH2D*> convertHistogram(TH2D* originalHist, const map<int, ChannelInfo>& channelMap, int rack) {
    vector<TH2D*> histograms;

    int nXBins = originalHist->GetNbinsX(); // Number of incidence angle bins
    int nYBins = originalHist->GetNbinsY(); // Number of channels

    // Iterate through X bins (one TH2D for each incidence angle)
    for (int xBin = 1; xBin <= nXBins; ++xBin) 
    {
        double angle = originalHist->GetXaxis()->GetBinCenter(xBin);  // Incidence angle
        string histName = Form("SlotCh_Angle_%.1f - %d", angle, rack);
        
        TH2D* slotChHist = new TH2D(histName.c_str(), 
                                    Form("Slot vs Ch (Angle %.1f);Slot;Ch", angle), 
                                    24, 0.5, 24+0.5,   // Example Slot range
                                    31, 0.5, 31+0.5); // Example Ch range

        // Iterate through Y bins (each corresponding to a channel number)
        for (int yBin = 1; yBin <= nYBins; ++yBin) 
        {
            int channel = originalHist->GetYaxis()->GetBinCenter(yBin); // Channel number
            double efficiency = originalHist->GetBinContent(xBin, yBin);

            // Check if channel exists in mapping
            if (channelMap.find(channel) != channelMap.end()) 
            {
                ChannelInfo info = channelMap.at(channel);
            	if (info.Rack == rack)
            	{
                	slotChHist->Fill(info.Slot, info.Ch, efficiency);
            	}
            }
        }

        histograms.push_back(slotChHist);
    }

    return histograms;
}


void drawFancyEffMap(TH2D* originalHist, int rack)
{

	std::string histName = originalHist->GetName();

	TCanvas* canvas = new TCanvas(histName.c_str(), histName.c_str(), 4096, 2160);
	canvas->SetGrid(0, 0);  // 0,0 means no grid in both X and Y axes

	TH2D *histData = (TH2D*)originalHist->Clone("histData");

	TH2D *histEmpty = (TH2D*)originalHist->Clone("histEmpty");
	histEmpty->Reset();  // Clears bin content while keeping axis and structure

	TH2D *histMissing = (TH2D*)originalHist->Clone("histMissing");
	histMissing->Reset();  // Clears bin content while keeping axis and structure

	// Fill the histograms
	for (int i = 1; i <= histData->GetNbinsX(); ++i) 
	{
	    for (int j = 1; j <= histData->GetNbinsY(); ++j) 
	    {
	        double value = histData->GetBinContent(i, j);

	        if (value == 0 && rack == 7)
	        {
	        	if ((j>=0 && j<=12) || (j>=16 && j<=28))
	        	{
	        		if (i == 2 || i == 5 || i == 8 || i== 11 || i ==19 || i == 22)
	        		{
	        			histMissing->SetBinContent(i, j, 1);
	        		}
	        		else histEmpty->SetBinContent(i, j, 1);
	        	}
	        	else histEmpty->SetBinContent(i, j, 1);
	        }

	       	if (value == 0 && rack == 8) 
	        {
				if ((j>=0 && j<=12) || (j>=16 && j<=28))
	        	{
	        		if (i == 5 || i == 8 || i== 11 || i ==14 || i == 17)
	        		{
	        			histMissing->SetBinContent(i, j, 1);
	        		}
	        		else histEmpty->SetBinContent(i, j, 1);
	        	}
	        	else histEmpty->SetBinContent(i, j, 1);	        
        	}

	    }
	}
	gPad->SetRightMargin(.15);	

	// Set up styles
	histEmpty->SetFillColor(kBlack);
	histEmpty->SetFillStyle(3004); // Cross-hatching pattern
	histEmpty->SetLineColor(kBlack); // Optional: border for clarity

	// Set up styles
	histMissing->SetFillColor(kRed);
	histMissing->SetFillStyle(3002); // Cross-hatching pattern
	histMissing->SetLineColor(kRed); // Optional: border for clarity

	// Draw histograms
	histData->GetXaxis()->SetTitle("Slot Number");
	histData->GetYaxis()->SetTitle("Channel Number");
	histData->GetZaxis()->SetTitle("Efficiency");

	histData->SetStats(kFALSE);  // Disable the stats box
	histData->GetXaxis()->SetTickLength(0);
	histData->GetYaxis()->SetTickLength(0);
	histData->Draw("COLZ");

	histEmpty->Draw("SAME BOX");  // Overlay empty bins with cross-hatching
	histMissing->Draw("SAME BOX");

	// Draw boxes around each bin
	for (int i = 1; i <= histData->GetNbinsX(); ++i) 
	{
	    for (int j = 1; j <= histData->GetNbinsY(); ++j) 
	    {
	        
	        // Bin edges
	        double xLow  = histData->GetXaxis()->GetBinLowEdge(i);
	        double xHigh = histData->GetXaxis()->GetBinUpEdge(i);
	        double yLow  = histData->GetYaxis()->GetBinLowEdge(j);
	        double yHigh = histData->GetYaxis()->GetBinUpEdge(j);

	        // Draw the box
	        TBox *box = new TBox(xLow, yLow, xHigh, yHigh);
	        box->SetLineColor(kBlack);   // Box outline color
	        box->SetLineWidth(1);        // Outline thickness
	        box->SetFillStyle(0);        // No fill (just the border)
	        box->Draw("same");
	    }
	}
	canvas->Update();

	canvas->Write(histName.c_str());
	string filename = histName + ".png";
	canvas->SaveAs(filename.c_str());
}


// void drawFitParameters(TF1 * fitFunction)
// {

// 	// 2. Extract fit parameters and errors
// 	double param1 = fitFunction->GetParameter(0); // Parameter 1 (amplitude)
// 	double param2 = fitFunction->GetParameter(1); // Parameter 2 (mean)
// 	double param3 = fitFunction->GetParameter(2); // Parameter 3 (sigma)

// 	double param1Error = fitFunction->GetParError(0); // Error for parameter 1
// 	double param2Error = fitFunction->GetParError(1); // Error for parameter 2
// 	double param3Error = fitFunction->GetParError(2); // Error for parameter 3



// 	// 3. Format the text for the label using #splitline{}

// 	string labelText1 =  "#splitline{Fit Results:}{Amplitude = " + Form("%.0f", param1) + "#pm" + param1Error;
// 	string labelText2 =  "#splitline{Mean = " + Form("%.0f", param2) + " #pm "+ param2Error +"}{Sigma = "+ param3 +" #pm "+param3Error +"}";
// 	string labelText  = "splitline{" + labelText1 +"}{"+labelText2+"}";


// 	// 4. Create a single TPaveLabel with the split text
// 	TPaveLabel *label = new TPaveLabel(0.6, 0.6, 0.9, 0.9, labelText.c_str(), "NDC");

// 	// 5. Customize the label (e.g., size, fill color)
// 	label->SetTextSize(0.04);
// 	label->SetFillColor(0);  // Transparent background
// 	label->SetBorderSize(1);

// 	// 6. Draw the label on the canvas
// 	label->Draw();
// }