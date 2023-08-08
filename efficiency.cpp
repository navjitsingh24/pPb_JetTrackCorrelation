#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TH3D.h"
#include "TVector.h"
#include "TDirectory.h"
#include "random_mixing.h"


void efficiency(TString input, TString output){

    float pi = TMath::Pi();

    TFile *input_file = new TFile(Form("%s",input.Data()));
    TTree *gen_tree = (TTree*)input_file->Get("HiGenParticleAna/hi");
    TTree *trk_tree = (TTree*)input_file->Get("ppTrack/trackTree");
    TTree *jet_tree = (TTree*)input_file->Get("ak4PFJetAnalyzer/t");
    TTree *pthat_tree = (TTree*)input_file->Get("hiEvtAnalyzer/HiTree");
    TTree *filter_tree = (TTree*)input_file->Get("skimanalysis/HltTree");
  
    // could also add tress instead of calling them individually inside the event loop
  
    //jet_tree->AddFriend(filter_tree);
    //jet_tree->AddFriend(trk_tree);
    //jet_tree->AddFriend(gen_tree);

    // defining the branches 

    int arr_size = 49999;

    float genjet_eta[arr_size];
    float genjet_pt[arr_size];
    float genjet_phi[arr_size];

    float recojet_eta[arr_size];
    float recojet_pt[arr_size];
    float recojet_phi[arr_size];

    float trk_eta[arr_size];
    float trk_pt[arr_size];
    float trk_pterr[arr_size];
    float trk_dxy[arr_size];
    float trk_dxyerr[arr_size];
    float trk_dz[arr_size];
    float trk_dzerr[arr_size];
    float trk_phi[arr_size];

    bool high_purity[arr_size];

    int  track_charge[arr_size];

    std::vector<float> *gen_pt = 0;
    std::vector<float> *gen_eta = 0;
    std::vector<float> *gen_phi = 0;
    std::vector<int> *gen_charge = 0;

    int ntrks, ngen, ngenjet, nrecojet;

    float pT_hat, evtweight, Vz;

    int pbeamscraping, ppaprimaryvertex, hbhenoiseresultrun2loose, phfcoinc, pvertexcutdz1p0;

    jet_tree->SetBranchAddress("ngen",&ngenjet);
    jet_tree->SetBranchAddress("nref",&nrecojet);

    pthat_tree->SetBranchAddress("pthat",&pT_hat);
    pthat_tree->SetBranchAddress("vz",&Vz);

    filter_tree->SetBranchAddress("pBeamScrapingFilter",&pbeamscraping);
    filter_tree->SetBranchAddress("pPAprimaryVertexFilter",&ppaprimaryvertex);
    filter_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&hbhenoiseresultrun2loose);
    filter_tree->SetBranchAddress("phfCoincFilter",&phfcoinc);
    filter_tree->SetBranchAddress("pVertexFilterCutdz1p0",&pvertexcutdz1p0);

    jet_tree->SetBranchAddress("genpt",&genjet_pt);
    jet_tree->SetBranchAddress("WTAgeneta",&genjet_eta);
    jet_tree->SetBranchAddress("WTAgenphi",&genjet_phi);
    jet_tree->SetBranchAddress("rawpt",&recojet_pt);
    jet_tree->SetBranchAddress("WTAeta",&recojet_eta);
    jet_tree->SetBranchAddress("WTAphi",&recojet_phi);

    trk_tree->SetBranchAddress("nTrk",&ntrks);
    trk_tree->SetBranchAddress("trkEta",&trk_eta);
    trk_tree->SetBranchAddress("trkPt",&trk_pt);
    trk_tree->SetBranchAddress("trkPtError",&trk_pterr);
    trk_tree->SetBranchAddress("trkDxy1",&trk_dxy);
    trk_tree->SetBranchAddress("trkDxyError1",&trk_dxyerr);
    trk_tree->SetBranchAddress("trkDz1",&trk_dz);
    trk_tree->SetBranchAddress("trkDzError1",&trk_dzerr);
    trk_tree->SetBranchAddress("trkPhi",&trk_phi);
    trk_tree->SetBranchAddress("trkCharge",&track_charge);
    trk_tree->SetBranchAddress("highPurity",&high_purity);

    gen_tree->SetBranchAddress("pt",&gen_pt);
    gen_tree->SetBranchAddress("eta",&gen_eta);
    gen_tree->SetBranchAddress("phi",&gen_phi);
    gen_tree->SetBranchAddress("chg",&gen_charge);

    int bins[3] ={12,16,51};
    double xmin[3]={-4.0,-pi,0.0};
    double xmax[3]={4.0,pi,500.0};

    int bins1[3] ={12,16,51};
    double xmin1[3]={-2.4,-pi,0.0};
    double xmax1[3]={2.4,pi,500.0};

    int npt_bins = 51;
    double pTbins[52] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 45.0, 60.0, 90.0, 120.0, 180.0, 300.0, 500.0};

    int neta_bins = 12;
    double etabins[13] = {-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0,2.4};

    int multiplicity_bins = 80;
    double multiplicity_step = 5.0;
    double multiplicitybins[81];
    for (int i=0;i<81;i++){
        multiplicitybins[i] = 0.0+i*multiplicity_step;
    }

    int dr_bins = 14;
    double drbins[15]={0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.0};

    int nphi_bins1 = 16;
    double phi_step1 = (2*pi)/16;
    double phibins1[17];
    for (int i = 0; i < 17; i++) {
        phibins1[i] = -pi + i*phi_step1;
    }


    THnSparseD *gen_trks = new THnSparseD("gen_trks","Gen;#eta;#phi;p_{T}",3,bins,xmin,xmax); 
    gen_trks->GetAxis(0)->Set(bins[0],etabins);
    gen_trks->GetAxis(1)->Set(bins[1],phibins1);
    gen_trks->GetAxis(2)->Set(bins[2],pTbins);
    gen_trks->Sumw2();
    THnSparseD *reco_trks = new THnSparseD("reco_trks","Reco;#eta;#phi;p_{T}",3,bins,xmin,xmax); 
    reco_trks->GetAxis(0)->Set(bins[0],etabins);
    reco_trks->GetAxis(1)->Set(bins[1],phibins1);
    reco_trks->GetAxis(2)->Set(bins[2],pTbins);
    reco_trks->Sumw2();
    THnSparseD *match_trks = new THnSparseD("match_trks","Match;#eta;#phi;p_{T}",3,bins,xmin,xmax); 
    match_trks->GetAxis(0)->Set(bins[0],etabins);
    match_trks->GetAxis(1)->Set(bins[1],phibins1);
    match_trks->GetAxis(2)->Set(bins[2],pTbins);
    match_trks->Sumw2();
    THnSparseD *fake_trks = new THnSparseD("fake_trks","Fake;#eta;#phi;p_{T}",3,bins,xmin,xmax); 
    fake_trks->GetAxis(0)->Set(bins[0],etabins);
    fake_trks->GetAxis(1)->Set(bins[1],phibins1);
    fake_trks->GetAxis(2)->Set(bins[2],pTbins);
    fake_trks->Sumw2();

    TH1D *pthat_before_w = new TH1D("pthat_before_w","Before Weighting;pthat",200,0.0,1000.0);
    pthat_before_w->Sumw2();
    TH1D *pthat_after_w = new TH1D("pthat_after_w","After Weighting;pthat",200,0.0,1000.0);
    pthat_after_w->Sumw2();

    TH1D *vz_dist = new TH1D("vz_dist",";V_{z}",60,-15.0,15.0);
    vz_dist->Sumw2();

    int nev = jet_tree->GetEntries();

    std::vector<float> dR;
    std::vector<int> min_gen_index;


    for (int i=0;i<nev;i++){ // event loop

        gen_tree->GetEntry(i);
        trk_tree->GetEntry(i);
        filter_tree->GetEntry(i);
        
        jet_tree->GetEntry(i);
        pthat_tree->GetEntry(i);

        if (pbeamscraping != 1) continue;
        if (ppaprimaryvertex != 1) continue;
        if (hbhenoiseresultrun2loose != 1) continue;
        if (phfcoinc != 1) continue;
        if (pvertexcutdz1p0 != 1) continue;

        if (Vz >15.0 || Vz < -15.0) continue;

        /*if(pT_hat > 15.0 && pT_hat <= 30.){evtweight = 1.0404701e-06 * 961104;}
		    else if(pT_hat > 30. && pT_hat <= 50.){evtweight = 7.7966624e-08 * 952110;}
		    else if(pT_hat > 50. && pT_hat <= 80.){evtweight = 1.0016052e-08 * 952554;}
		    else if(pT_hat > 80. && pT_hat <= 120.){evtweight = 1.3018269e-09 * 996844;}
		    else if(pT_hat > 120. && pT_hat <= 170.){evtweight = 2.2648493e-10 * 964681;}
		    else if(pT_hat > 170. && pT_hat <= 220.){evtweight = 4.0879112e-11 * 999260;}
		    else if(pT_hat > 220. && pT_hat <= 280.){evtweight = 1.1898939e-11 * 964336;}
		    else if(pT_hat > 280. && pT_hat <= 370.){evtweight = 3.3364433e-12 * 995036;}
		    else if(pT_hat > 370. && pT_hat <= 460.){evtweight = 7.6612402e-13 * 958160;}
		    else if(pT_hat > 460. && pT_hat <= 540.){evtweight = 2.1341026e-13 * 981427;}
		    else if(pT_hat > 540.){evtweight = 7.9191586e-14 * 1000000;}
        evtweight = (float) evtweight / nev;*/

        evtweight = 1.0;

        pthat_before_w->Fill(pT_hat);
        pthat_after_w->Fill(pT_hat,evtweight);

        vz_dist->Fill(Vz,evtweight);

        int multiplicity = get_Ntrkoff(ntrks, trk_eta, trk_pt, track_charge, high_purity, trk_pterr, trk_dxy, trk_dxyerr, trk_dz, trk_dzerr);
        
    	  std::vector<TVector3> matchreco;
    	  std::vector<TVector3> gentrk;
        std::vector<TVector3> faketrk;
    	  std::vector<TVector3> recotrk;
        std::vector<TVector3> matchgen;

        ngen = gen_eta->size();

        for (int j=0;j<ngenjet;j++){ // jet loop
            //h_recojet_spectra->Fill(recojet_eta[j], recojet_phi[j], recojet_pt[j],evtweight);
            if (genjet_pt[j]<=60.0) continue;
            if (fabs(genjet_eta[j])>1.6) continue;
            TVector3 JetTracks;
	 		      JetTracks.SetPtEtaPhi(genjet_pt[j], genjet_eta[j], genjet_phi[j]);
	 		      jettrk.push_back(JetTracks);
            jet_weight.push_back(1.0);
            jetpt_histo->Fill(genjet_pt[j],evtweight);
        }

        for (int q=0;q<ntrks;q++){//reco loop
            if (track_charge[q]==0) continue;
            if (high_purity[q]!=1) continue;
            if (fabs(trk_dxy[q]/trk_dxyerr[q])>=3.0) continue;
            if (fabs(trk_dz[q]/trk_dzerr[q])>=3.0) continue;
            if (fabs(trk_pterr[q]/trk_pt[q])>=0.1) continue;
            if (trk_pt[q]<=0.4) continue;
            if (fabs(trk_eta[q])>2.4) continue;
            TVector3 RecoTracks;
            RecoTracks.SetPtEtaPhi(trk_pt[q],trk_eta[q],trk_phi[q]);
            recotrk.push_back(RecoTracks);
            float trkweight;
            if (trk_pt[q]<=8.0){
               trkweight = getTrkCorrWeight(hijing_effcy,hijing_fkrate,trk_pt[q],trk_eta[q],trk_phi[q]);
            }
            else if (trk_pt[q]>8.0){
                trkweight = getTrkCorrWeight(pythia_effcy,pythia_fkrate,trk_pt[q],trk_eta[q],trk_phi[q]);
            }
            trk_weight.push_back(trkweight);
            trk_weight_one.push_back(1.0);
        } 

        for (int x=0;x<ngen;x++){ //gen loop 
            if (gen_charge->at(x)==0) continue;
            if (gen_pt->at(x) <= 0.4) continue;
            if (fabs(gen_eta->at(x))>2.4) continue;
            TVector3 GenTracks;
	 		     GenTracks.SetPtEtaPhi(gen_pt->at(x), gen_eta->at(x), gen_phi->at(x));
	 	       gentrk.push_back(GenTracks);
		    }

        //if (recotrk.size() ==0 || gentrk.size()==0 ) continue;

        for (int p=0;p<gentrk.size(); p++){
            double gentrk_variables[3]= {gentrk[p].Eta(),gentrk[p].Phi(),gentrk[p].Pt()};
            gen_trks->Fill(gentrk_variables,evtweight);
        }
        for (int w=0;w<recotrk.size();w++){
            double recotrk_variables[3] = {recotrk[w].Eta(),recotrk[w].Phi(),recotrk[w].Pt()};
            reco_trks_eff_corr->Fill(recotrk_variables,(evtweight*trk_weight.at(w))); 
        }
      
        for (int j=0;j<ntrks;j++){ //reco loop
            dR.clear();
            min_gen_index.clear();
            if (track_charge[j]==0) continue;
            if (high_purity[j]!=1) continue;
            if (fabs(trk_dxy[j]/trk_dxyerr[j])>=3.0) continue;
            if (fabs(trk_dz[j]/trk_dzerr[j])>=3.0) continue;
            if (fabs(trk_pterr[j]/trk_pt[j])>=0.1) continue;
            if (trk_pt[j]<=0.4) continue;
            if (fabs(trk_eta[j])>2.4) continue;
            for (int k=0;k<ngen;k++){ // gen loop
                if (gen_charge->at(k)==0) continue;
                if (gen_pt->at(k)<= 0.4) continue;
                if (fabs(gen_eta->at(k))>2.4) continue;
                float deltaeta = trk_eta[j]-gen_eta->at(k);
                float deltaphi = trk_phi[j]-gen_phi->at(k);
                if (deltaphi > pi){deltaphi+= -2*pi;}
                else if (deltaphi < -pi){deltaphi+= 2*pi;}
                float deltarsq = (deltaeta*deltaeta) + (deltaphi*deltaphi);
                float deltaR = TMath::Sqrt(deltarsq);
                dR.push_back(deltaR);
                min_gen_index.push_back(k);
            }

            //if (dR.size()==0) continue; 

            float drcut = *min_element(dR.begin(),dR.end());
            int drcut_elmnt = min_element(dR.begin(),dR.end()) -dR.begin();
            int dr_min_index = min_gen_index.at(drcut_elmnt);
            float dptcut = fabs(trk_pt[j]-gen_pt->at(dr_min_index))/gen_pt->at(dr_min_index);
            TVector3 MatchedRecoTracks;
	 		      MatchedRecoTracks.SetPtEtaPhi(trk_pt[j], trk_eta[j], trk_phi[j]);
            TVector3 MatchedGenTracks;
	 		      MatchedGenTracks.SetPtEtaPhi(gen_pt->at(dr_min_index), gen_eta->at(dr_min_index), gen_phi->at(dr_min_index));
            // matching condition
            if (drcut < 0.1 && dptcut < 0.1){
                matchreco.push_back(MatchedRecoTracks);
                matchgen.push_back(MatchedGenTracks);
            } // fake trks
            else if (drcut >= 0.1 || dptcut >= 0.1 ){
                faketrk.push_back(MatchedRecoTracks);
            }
        } 

      for (int p=0;p<matchreco.size(); p++){
          double matchtrk_variables[3]= {matchreco[p].Eta(),matchreco[p].Phi(),matchreco[p].Pt()};
          match_trks->Fill(matchtrk_variables,evtweight);
      }
      for (int w=0;w<faketrk.size();w++){
          double faketrk_variables[3] = {faketrk[w].Eta(),faketrk[w].Phi(),faketrk[w].Pt()};
          fake_trks->Fill(faketrk_variables,evtweight); 
      }
      
    
    }

    TFile *out_file = new TFile(Form("%s",output.Data()),"recreate");
    out_file->cd();
    gen_trks->Write();
    reco_trks->Write();
    match_trks->Write();
    fake_trks->Write();
    pthat_before_w->Write();
    pthat_after_w->Write();
    vz_dist->Write();

}
