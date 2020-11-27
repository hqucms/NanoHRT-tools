#include "setTDRStyle.h"
#include <sstream>
#include <fstream>
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include <vector>
#include "Rtypes.h"
#include "TColor.h"
#include "TVectorF.h"
#include <cstdlib>

void setTDRStyle();
TH1F *create1Dhisto(TString sample,TTree *tree,TString intLumi,TString cuts,TString branch,int bins,float xmin,float xmax,
                    bool useLog,int color, int style,TString name,bool norm,bool data);
TCanvas *c_comp1Dhistos(TString c_name, TString lumi, std::vector<TString> names, std::vector<TTree*> trees, std::vector<bool> isdata, std::vector<TString> vars, std::vector<TString> cuts, 
			std::vector<int> colors, std::vector<int> styles, TString xname, TString yname, 
		        int nbins, float xmin, float xmax, bool uselog, bool norm);
TCanvas *c_comp1Dhistos(TString c_name, TString lumi, std::vector<TString> names, std::vector<TTree*> trees, std::vector<bool> isdata, TString var, std::vector<TString> cuts, 
			std::vector<int> colors, std::vector<int> styles, TString xname, TString yname, 
		        int nbins, float xmin, float xmax, bool uselog, bool norm);
TCanvas *c_ratioOf1Dhistos(TString c_name, TString lumi, std::vector<TString> namest, std::vector<TTree*> trees, std::vector<bool> isdata, TString var, 
			   std::vector<TString> namesc, std::vector<TString> cuts,
                           TString xname, TString yname, int nbins, float xmin, float xmax, bool uselog, bool norm);

TH1D *hrwgt = nullptr;

double rewgtfunc(double x){
  if (!hrwgt) return 1;
  return hrwgt->GetBinContent(hrwgt->FindFixBin(x));
}

void compDiscInProxyVsSignalJets(TString bkgSample, TString year, TString cat, TString score) {

  TH1::SetDefaultSumw2(kTRUE);
  setTDRStyle();
  gROOT->SetBatch(false);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(1);

  float intLumi;
  if (year == "2016") { intLumi = 36.8; }
  ostringstream tmpLumi;
  tmpLumi << intLumi;
  TString lumi = tmpLumi.str();


  TFile *f_signal = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200316/ggH-bb_tree.root" , "READONLY");
  TFile *f_signal_b;
  if (cat == "bb") { f_signal_b = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200330/vhtobb_tree.root" , "READONLY"); }
  if (cat == "cc") { f_signal_b = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200330/vhtocc_tree.root" , "READONLY"); }
  //  TFile *f_proxy  = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200312/qcd-mg_tree.root" , "READONLY"); // sm_tree.root dytonunu_tree.root
  //TFile *f_proxy  = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200316/znunu_tree.root" , "READONLY");
  TFile *f_proxy  = TFile::Open("/eos/uscms/store/user/lpcjme/noreplica/loukas/particlenet/trees/inclusive/2016/20200330/dytonunu_tree.root" , "READONLY");

  TTree *t_signal   = (TTree*)f_signal->Get("Events");
  TTree *t_signal_b = (TTree*)f_signal_b->Get("Events");  
  TTree *t_proxy    = (TTree*)f_proxy->Get("Events");
  
  std::vector<TTree*> trees;  trees.size();
  std::vector<bool>   isdata; isdata.size(); 
  std::vector<int>    colors; colors.size(); 
  std::vector<int>    styles; styles.size();
  //trees.push_back(t_signal);   isdata.push_back(false); colors.push_back(1); styles.push_back(1);
  trees.push_back(t_signal_b); isdata.push_back(false); colors.push_back(1); styles.push_back(1);
  trees.push_back(t_signal_b); isdata.push_back(false); colors.push_back(1); styles.push_back(2);
  trees.push_back(t_proxy);    isdata.push_back(false); colors.push_back(2); styles.push_back(1);
  trees.push_back(t_proxy);    isdata.push_back(false); colors.push_back(2); styles.push_back(2);
  trees.push_back(t_proxy);    isdata.push_back(false); colors.push_back(4); styles.push_back(1);
  trees.push_back(t_proxy);    isdata.push_back(false); colors.push_back(4); styles.push_back(2);
  trees.push_back(t_proxy);    isdata.push_back(false); colors.push_back(6); styles.push_back(1);
  trees.push_back(t_proxy);    isdata.push_back(false); colors.push_back(6); styles.push_back(2);

  std::vector<TString> names; names.clear();
  //names.push_back("signal-ggh");
  names.push_back("signal_vh");
  names.push_back("signal-b-vh");
  names.push_back("znunu-a");
  names.push_back("znunu-rwgt-a");
  names.push_back("znunu-b");
  names.push_back("znunu-rwgt-b");
  names.push_back("znunu-c");
  names.push_back("znunu-rwgt-c");

  TString c_match; 
  if (cat == "bb") { c_match = "fj_1_nbhadrons>=2"; } 
  if (cat == "cc") { c_match = "fj_1_nbhadrons==0 && fj_1_nchadrons>=2"; }
  
  TString c_inc     = "fj_1_pt>300. && fj_1_pt<400. && abs(fj_1_eta)<2.4 && fj_1_sdmass>50. && fj_1_sdmass<200. && ht>200. && fj_1_tau21<1.1 && (fj_1_sdmass/fj_1_energy)>0.";
  TString c_fjsv2   = "fj_1_nsv>1";
  TString c_sjsv2   = "fj_1_sj1_nsv>0 && fj_1_sj2_nsv>0";
  TString c_sjsv2_t = "(fj_1_sj1_nsv>0 && fj_1_sj1_sv1_ntracks>2 && abs(fj_1_sj1_sv1_dxy)<3 && fj_1_sj1_sv1_dlensig>4.) && (fj_1_sj2_nsv>0 && fj_1_sj2_sv1_ntracks>2 && abs(fj_1_sj2_sv1_dxy)<3 && fj_1_sj2_sv1_dlensig>4.)";

  //TString c_incl         = c_inc+" && "+c_match+" && "+c_sjsv2_t;
  TString c_incl         = c_inc+" && "+c_match;
  TString c_incl_fjsv2   = c_inc+" && "+c_fjsv2+" && "+c_match;
  TString c_incl_sjsv2   = c_inc+" && "+c_sjsv2+" && "+c_match;
  TString c_incl_sjsv2_t = c_inc+" && "+c_sjsv2_t+" && "+c_match;

  std::vector<TString> cuts; cuts.clear();
  cuts.push_back(c_incl+" && fj_1_isH<0.6");
  cuts.push_back(c_incl_sjsv2_t+" && fj_1_isH<0.6");
  cuts.push_back(c_incl+" && "+c_match);
  cuts.push_back(c_incl+" && "+c_match);
  cuts.push_back(c_incl_sjsv2+" && "+c_match);
  cuts.push_back(c_incl_sjsv2+" && "+c_match);
  cuts.push_back(c_incl_sjsv2_t+" && "+c_match);
  cuts.push_back(c_incl_sjsv2_t+" && "+c_match);

  // reweighting of two distributions
  TH1F *h_test_1 = create1Dhisto("signal",t_signal_b,lumi,c_incl+" && "+c_match+" && fj_1_isH<0.6","fj_1_tau21",20,0.,1.,false,1,1,"h_test_1",true,false); h_test_1->SetFillColor(0);
  TH1F *h_test_2 = create1Dhisto("signal",t_proxy,lumi,c_incl+" && "+c_match,"fj_1_tau21",20,0.,1.,false,2,2,"h_test_2",true,false); h_test_2->SetFillColor(0);

  //  TH1F *h_test_1 = create1Dhisto("signal",t_signal_b,lumi,c_incl+" && "+c_match+" && fj_1_isH<0.8","fj_1_sdmass",20,50.,250.,false,1,1,"h_test_1",true,false); h_test_1->SetFillColor(0);
  //  TH1F *h_test_2 = create1Dhisto("signal",t_proxy,lumi,c_incl+" && "+c_match,"fj_1_sdmass",20,50.,250.,false,2,2,"h_test_2",true,false); h_test_2->SetFillColor(0);
  
  //TH1F *h_test_1 = create1Dhisto("signal",t_signal_b,lumi,c_incl+" && fj_1_isH<0.8","sqrt((fj_1_sj1_eta-fj_1_sj2_eta)*(fj_1_sj1_eta-fj_1_sj2_eta) + (fj_1_sj1_eta-fj_1_sj2_phi)*(fj_1_sj1_eta-fj_1_sj2_phi))",35,0.,7.,false,1,1,"h_test_1",true,false); h_test_1->SetFillColor(0);
  //TH1F *h_test_2 = create1Dhisto("signal",t_proxy,lumi,c_incl+" && "+c_match,"sqrt((fj_1_sj1_eta-fj_1_sj2_eta)*(fj_1_sj1_eta-fj_1_sj2_eta) + (fj_1_sj1_eta-fj_1_sj2_phi)*(fj_1_sj1_eta-fj_1_sj2_phi))",35,0.,7.,false,2,2,"h_test_2",true,false); h_test_2->SetFillColor(0);
  hrwgt = (TH1D*)h_test_1->Clone("hrwgt"); hrwgt->Divide(h_test_1,h_test_2);
  /*
  if (cat == "bb") {  
    TString var_fj_1_ParticleNetMD_XbbVsQCD ="fj_1_ParticleNetMD_XbbVsQCD";
    TCanvas *c_fj_1_ParticleNetMD_XbbVsQCD = c_comp1Dhistos("fj_1_ParticleNetMD_XbbVsQCD_"+cat,lumi,names,trees,isdata,var_fj_1_ParticleNetMD_XbbVsQCD,cuts,colors,styles,
							    "ParticleNetMD_XbbVsQCD","a. u.",20,0.,1.,false,true);
  }

  if (cat == "cc") {
  TString var_fj_1_ParticleNetMD_XccVsQCD ="fj_1_ParticleNetMD_XccVsQCD";
  TCanvas *c_fj_1_ParticleNetMD_XccVsQCD = c_comp1Dhistos("fj_1_ParticleNetMD_XccVsQCD_"+cat,lumi,names,trees,isdata,var_fj_1_ParticleNetMD_XccVsQCD,cuts,colors,styles,
							  "ParticleNetMD_XccVsQCD","a. u.",20,0.,1.,false,true);
  }
  */

  /*
  TString var_fj_1_ParticleNetMD_bbVsLight ="fj_1_ParticleNetMD_bbVsLight";
  TCanvas *c_fj_1_ParticleNetMD_bbVsLight = c_comp1Dhistos("fj_1_ParticleNetMD_bbVsLight_"+cat,lumi,names,trees,isdata,var_fj_1_ParticleNetMD_bbVsLight,cuts,colors,styles,
							  "ParticleNetMD_bbVsLight","a. u.",20,0.,1.,false,true);

  TString var_fj_1_ParticleNetMD_ccVsLight ="fj_1_ParticleNetMD_ccVsLight";
  TCanvas *c_fj_1_ParticleNetMD_ccVsLight = c_comp1Dhistos("fj_1_ParticleNetMD_ccVsLight_"+cat,lumi,names,trees,isdata,var_fj_1_ParticleNetMD_ccVsLight,cuts,colors,styles,
							  "ParticleNetMD_ccVsLight","a. u.",20,0.,1.,false,true);

  TString var_fj_1_DeepAK8MD_ZHbbvsQCD ="fj_1_DeepAK8MD_ZHbbvsQCD";
  TCanvas *c_fj_1_DeepAK8MD_ZHbbvsQCD = c_comp1Dhistos("fj_1_DeepAK8MD_ZHbbvsQCD_"+cat,lumi,names,trees,isdata,var_fj_1_DeepAK8MD_ZHbbvsQCD,cuts,colors,styles,
						       "DeepAK8MD_ZHbbvsQCD","a. u.",20,0.,1.,false,true);

  TString var_fj_1_DeepAK8MD_ZHccvsQCD ="fj_1_DeepAK8MD_ZHccvsQCD";
  TCanvas *c_fj_1_DeepAK8MD_ZHccvsQCD = c_comp1Dhistos("fj_1_DeepAK8MD_ZHccvsQCD_"+cat,lumi,names,trees,isdata,var_fj_1_DeepAK8MD_ZHccvsQCD,cuts,colors,styles,
						       "DeepAK8MD_ZHccvsQCD","a. u.",20,0.,1.,false,true);

  TString var_fj_sj1_pt ="fj_1_sj1_pt";
  TCanvas *c_fj_sj1_pt = c_comp1Dhistos("fj_1_sj1_pt_"+cat,lumi,names,trees,isdata,var_fj_sj1_pt,cuts,colors,styles,"p_{T}(sj1) [GeV]","a. u.",40,0.,800.,false,true);
  
  TString var_fj_sj2_pt ="fj_1_sj2_pt";
  TCanvas *c_fj_sj2_pt = c_comp1Dhistos("fj_1_sj2_pt_"+cat,lumi,names,trees,isdata,var_fj_sj2_pt,cuts,colors,styles,"p_{T}(sj2) [GeV]","a. u.",40,0.,800.,false,true);
 
  TString var_fj_sj2ovsj1_pt ="fj_1_sj2_pt/fj_1_sj1_pt";
  TCanvas *c_fj_sj2ovsj1_pt = c_comp1Dhistos("fj_1_sj2ovsj1_pt_"+cat,lumi,names,trees,isdata,var_fj_sj2ovsj1_pt,cuts,colors,styles,"p_{T}(sj2)/p_{T}(sj1)","a. u.",24,0.,1.2,false,true);
  
  TString var_fj_tau21 ="fj_1_tau21";
  TCanvas *c_fj_tau21 = c_comp1Dhistos("fj_1_tau21_"+cat,lumi,names,trees,isdata,var_fj_tau21,cuts,colors,styles,"#tau_{21}","a. u.",10,0.,1.,false,true);
  
  TString var_fj_nsv ="fj_1_nsv";
  TCanvas *c_fj_nsv = c_comp1Dhistos("fj_1_nsv_"+cat,lumi,names,trees,isdata,var_fj_nsv,cuts,colors,styles,"N_{SV}(fj)","a. u.",7,-0.5,6.5,false,true);

  TString var_fj_ntracks ="fj_1_ntracks";
  TCanvas *c_fj_ntracks = c_comp1Dhistos("fj_1_ntracks_"+cat,lumi,names,trees,isdata,var_fj_ntracks,cuts,colors,styles,"N_{tracks}(fj)","a. u.",19,-0.5,18.5,false,true);

  TString var_fj_btagcsvv2 ="fj_1_btagcsvv2";
  TCanvas *c_fj_btagcsvv2 = c_comp1Dhistos("fj_1_btagcsvv2_"+cat,lumi,names,trees,isdata,var_fj_btagcsvv2,cuts,colors,styles,"CSVv2(fj)","a. u.",10,0.,1.,false,true);
  
  TString var_ht ="ht";
  TCanvas *c_ht = c_comp1Dhistos("ht_"+cat,lumi,names,trees,isdata,var_ht,cuts,colors,styles,"ht","a. u.",50,0.,2000.,false,true);

  //  TString var_nlep ="nlep";
  //  TCanvas *c_nlep = c_comp1Dhistos("nlep_"+cat,lumi,names,trees,isdata,var_nlep,cuts,colors,styles,"nlep","a. u.",6,-0.5,5.5,false,true);
  
  TString var_sj12dr ="sqrt((fj_1_sj1_eta-fj_1_sj2_eta)*(fj_1_sj1_eta-fj_1_sj2_eta) + (fj_1_sj1_eta-fj_1_sj2_phi)*(fj_1_sj1_eta-fj_1_sj2_phi))";
  TCanvas *c_sj12dr = c_comp1Dhistos("fj_sj12dr_"+cat,lumi,names,trees,isdata,var_sj12dr,cuts,colors,styles,"#DeltaR(sj12)","a. u.",35,0.,7.,false,true);
  */
  /*
  TString var_fj_ntracks ="fj_1_ntracks";
  TCanvas *c_fj_ntracks = c_comp1Dhistos("fj_1_ntracks_"+cat,lumi,names,trees,isdata,var_fj_ntracks,cuts,colors,styles,"N_{tracks} (fj)","a. u.",31,-0.5,30.5,false,true);

  TString var_fj_sdovrawmass ="fj_1_sdmass/fj_1_rawmass";
  TCanvas *c_fj_sdovrawmass = c_comp1Dhistos("fj_1_sdovrawmass_"+cat,lumi,names,trees,isdata,var_fj_sdovrawmass,cuts,colors,styles,"m_{sd}/m_{raw} (fj)","a. u.",10,-0.,1.,false,true);

  TString var_fj_rawmassovenergy ="fj_1_rawmass/fj_1_energy";
  TCanvas *c_fj_rawmassovenergy = c_comp1Dhistos("fj_1_rawmassovenergy_"+cat,lumi,names,trees,isdata,var_fj_rawmassovenergy,cuts,colors,styles,"m_{raw}/E (fj)","a. u.",10,-0.,1.,false,true);

  TString var_fj_sdmassovenergy ="fj_1_sdmass/fj_1_energy";
  TCanvas *c_fj_sdmassovenergy = c_comp1Dhistos("fj_1_sdmassovenergy_"+cat,lumi,names,trees,isdata,var_fj_sdmassovenergy,cuts,colors,styles,"m_{sd}/E (fj)","a. u.",10,-0.,1.,false,true);
  */


  {
    std::vector<TString> namest;  namest.push_back("znunu"); namest.push_back("signal");
    std::vector<TString> namesc;  namesc.push_back("300to400"); namesc.push_back("400to500"); namesc.push_back("500toInf");
    std::vector<TTree*>  treesr;  treesr.push_back(t_proxy); treesr.push_back(t_signal_b);
    std::vector<bool>   isdatar; isdatar.push_back(false); isdatar.push_back(false);
    std::vector<TString>  cutsr;   
    cutsr.push_back(c_incl+"&& fj_1_pt>300. && fj_1_pt<400."); cutsr.push_back(c_incl+"&& fj_1_pt>400. && fj_1_pt<500."); cutsr.push_back(c_incl+"&& fj_1_pt>500. && fj_1_pt<50000.");
    TCanvas *c_ratio_fj_1_ntracks =c_ratioOf1Dhistos("r_fj_1_ntracks_"+cat,lumi,namest,treesr,isdatar,"fj_1_ntracks",namesc,cutsr,"N_{trks} (fj)","Ratio",31,-0.5,30.5,false,true);
  }

  
  /*
  TH1F *h_signal_incl    = create1Dhisto("signal",t_signal,lumi,c_incl+" && fj_1_isH<0.8",score,50,0.,1.,false,1,1,"h_signal_incl",true,false);            h_signal_incl->SetFillColor(0);
  // TH1F *h_signal_fjsv2   = create1Dhisto("signal",t_signal,lumi,c_incl_fjsv2+" && fj_1_isH<0.8",score,50,0.,1.,false,2,1,"h_signal_fjsv2",true,false);     h_signal_fjsv2->SetFillColor(0);
  //TH1F *h_signal_sjsv2   = create1Dhisto("signal",t_signal,lumi,c_incl_sjsv2+" && fj_1_isH<0.8",score,50,0.,1.,false,4,1,"h_signal_sjsv2",true,false);     h_signal_sjsv2->SetFillColor(0);
  TH1F *h_signal_sjsv2   = create1Dhisto("signal",t_signal,lumi,c_incl_sjsv2+" && fj_1_isH<0.8 && fj_1_sj1_nbhadrons>=1 && fj_1_sj1_nchadrons>=1 && fj_1_sj2_nbhadrons>=1 && fj_1_sj2_nchadrons>=1",score,50,0.,1.,false,4,1,"h_signal_sjsv2",true,false);     h_signal_sjsv2->SetFillColor(0);
  TH1F *h_signal_sjsv2_t = create1Dhisto("signal",t_signal,lumi,c_incl_sjsv2_t+" && fj_1_isH<0.8",score,50,0.,1.,false,6,1,"h_signal_sjsv2_t",true,false); h_signal_sjsv2_t->SetFillColor(0);
  

  // fj_2_ParticleNetMD_XbbVsQCD
  TH1F *h_proxy_incl    = create1Dhisto("qcd",t_proxy,lumi,c_incl+" && fj_1_partonflavour==21",score,50,0.,1.,false,1,2,"h_proxy_incl",true,false);            h_proxy_incl->SetFillColor(0);
  //  TH1F *h_proxy_fjsv2   = create1Dhisto("qcd",t_proxy,lumi,c_inc+" && fj_1_nbhadrons==0",score,50,0.,1.,false,2,2,"h_proxy_fjsv2",true,false);     h_proxy_fjsv2->SetFillColor(0);
  //TH1F *h_proxy_sjsv2   = create1Dhisto("qcd",t_proxy,lumi,c_incl_sjsv2,score,50,0.,1.,false,4,2,"h_proxy_sjsv2",true,false);     h_proxy_sjsv2->SetFillColor(0);
  TH1F *h_proxy_sjsv2   = create1Dhisto("qcd",t_proxy,lumi,c_incl+" && fj_1_partonflavour==21 &&  fj_1_sj1_nbhadrons>=1 && fj_1_sj1_nchadrons>=1 && fj_1_sj2_nbhadrons>=1 && fj_1_sj2_nchadrons>=1",score,50,0.,1.,false,4,2,"h_proxy_sjsv2",true,false);     h_proxy_sjsv2->SetFillColor(0);   
  TH1F *h_proxy_sjsv2_t = create1Dhisto("qcd",t_proxy,lumi,c_incl_sjsv2_t+" && fj_1_partonflavour==21",score,50,0.,1.,false,6,2,"h_proxy_sjsv2_t",true,false); h_proxy_sjsv2_t->SetFillColor(0);
*/
  
  /*
  TH1F *h_proxy_incl    = create1Dhisto("qcd",t_proxy,lumi,c_inc_b+" && fj_2_nbhadrons>=2","fj_1_ParticleNetMD_XbbVsQCD",50,0.,1.,false,1,2,"h_proxy_incl",true,false);            h_proxy_incl->SetFillColor(0);
  TH1F *h_proxy_incl_b    = create1Dhisto("qcd",t_proxy,lumi,c_inc_b+" && fj_2_nbhadrons==0 && fj_2_nchadrons==0","fj_1_ParticleNetMD_XbbVsQCD",50,0.,1.,false,2,2,"h_proxy_incl_b",true,false);            h_proxy_incl_b->SetFillColor(0);*/

  /*
  TCanvas *c_comp = new TCanvas("c_comp","c_comp",500,500);
  h_signal_incl->GetXaxis()->SetTitle(score);
  h_signal_incl->GetYaxis()->SetTitle("a. u.");
  h_signal_incl->GetYaxis()->SetRangeUser(0.,1.);
  h_signal_incl->Draw("HIST E0");
  //h_signal_fjsv2->Draw("HIST E0 sames");
  h_signal_sjsv2->Draw("HIST E0 sames");
  h_signal_sjsv2_t->Draw("HIST E0 sames");
  h_proxy_incl->Draw("HIST E0 sames");
  //h_proxy_fjsv2->Draw("HIST E0 sames");
  h_proxy_sjsv2->Draw("HIST E0 sames");
  h_proxy_sjsv2_t->Draw("HIST E0 sames");
  c_comp->RedrawAxis();
  */
  
  //TH1F *h_test_1 = create1Dhisto("qcd",t_signal,lumi,c_incl+" && fj_1_isH<0.8","fj_1_sj1_pt",100,0.,1000.,false,1,1,"h_test_1",true,false); h_test_1->SetFillColor(0);
  //TH1F *h_test_2 = create1Dhisto("qcd",t_proxy,lumi,c_incl+" && fj_1_nbhadrons>=2","fj_1_sj1_pt",100,0.,1000.,false,2,2,"h_test_2",true,false); h_test_2->SetFillColor(0);
  //  TH1F *h_test_1 = create1Dhisto("qcd",t_signal,lumi,c_incl+" && fj_1_isH<0.8","fj_1_sj2_pt/fj_1_sj1_pt",100,0.,1.,false,1,1,"h_test_1",true,false); h_test_1->SetFillColor(0);
  //  TH1F *h_test_2 = create1Dhisto("qcd",t_proxy,lumi,c_incl+" && fj_1_nbhadrons>=2","fj_1_sj2_pt/fj_1_sj1_pt",100,0.,1.,false,2,2,"h_test_2",true,false); h_test_2->SetFillColor(0);


}



TCanvas *c_comp1Dhistos(TString c_name, TString lumi, std::vector<TString> names, std::vector<TTree*> trees, std::vector<bool> isdata, std::vector<TString> vars, std::vector<TString> cuts, 
			std::vector<int> colors, std::vector<int> styles, TString xname, TString yname, 
		        int nbins, float xmin, float xmax, bool uselog, bool norm) {

  std::vector<TH1F*> histos; histos.clear(); 
  for (int itree=0; itree<trees.size(); ++itree) {
    TH1F *h_tmp = create1Dhisto(names[itree],trees[itree],lumi,cuts[itree],vars[itree],nbins,xmin,xmax,uselog,colors[itree],styles[itree],"h_"+names[itree]+"_"+c_name,norm,isdata[itree]); 
    h_tmp->SetFillColor(0);
    histos.push_back(h_tmp);
  }

  TCanvas *c_comp1Dhistos = new TCanvas("c_comp1Dhistos_"+c_name,"c_comp1Dhistos_"+c_name,500,500);
  for (int ihisto=0; ihisto<histos.size(); ++ihisto) {
    if (ihisto==0) {
      histos[ihisto]->GetXaxis()->SetTitle(xname);
      histos[ihisto]->GetYaxis()->SetTitle(yname);
      histos[ihisto]->GetYaxis()->SetRangeUser(0.,1.);
      histos[ihisto]->Draw("HIST E0");
    }
    else {
      histos[ihisto]->Draw("HIST E0 sames");
    }
  }

  return c_comp1Dhistos;
}


TCanvas *c_comp1Dhistos(TString c_name, TString lumi, std::vector<TString> names, std::vector<TTree*> trees, std::vector<bool> isdata, TString var, std::vector<TString> cuts, 
			std::vector<int> colors, std::vector<int> styles, TString xname, TString yname, 
		        int nbins, float xmin, float xmax, bool uselog, bool norm) {

  std::vector<TH1F*> histos; histos.clear(); 
  for (int itree=0; itree<trees.size(); ++itree) {
    TH1F *h_tmp = create1Dhisto(names[itree],trees[itree],lumi,cuts[itree],var,nbins,xmin,xmax,uselog,colors[itree],styles[itree],"h_"+names[itree]+"_"+c_name,norm,isdata[itree]); 
    h_tmp->SetFillColor(0);
    histos.push_back(h_tmp);
  }
  
  // get histo max value
  float maxYval = 0.;
  for (int ihisto=0; ihisto<histos.size(); ++ihisto) { 
    float maxYval_tmp = histos[ihisto]->GetBinContent(histos[ihisto]->GetMaximumBin());
    if (maxYval_tmp>maxYval) { maxYval = maxYval_tmp; }
  }

  TCanvas *c_comp1Dhistos = new TCanvas("c_comp1Dhistos_"+c_name,"c_comp1Dhistos_"+c_name,500,500);
  for (int ihisto=0; ihisto<histos.size(); ++ihisto) {
    if (ihisto==0) {
      histos[ihisto]->GetXaxis()->SetTitle(xname);
      histos[ihisto]->GetYaxis()->SetTitle(yname);
      histos[ihisto]->GetYaxis()->SetRangeUser(0.,1.3*maxYval);
      histos[ihisto]->Draw("HIST E0");
    }
    else {
      histos[ihisto]->Draw("HIST E0 sames");
    }
  }
  //gSystem->cd("/uscms_data/d3/loukas/ParticleNetInData/CMSSW_10_2_18/src/PhysicsTools/NanoHRTTools/plotting/plots/particlenet/");
  c_comp1Dhistos->Print("/uscms_data/d3/loukas/ParticleNetInData/CMSSW_10_2_18/src/PhysicsTools/NanoHRTTools/plotting/plots/particlenet/c_comp1Dhistos_"+c_name+".pdf");
  return c_comp1Dhistos;
}

TCanvas *c_ratioOf1Dhistos(TString c_name, TString lumi, std::vector<TString> namest, std::vector<TTree*> trees, std::vector<bool> isdata, TString var, 
			   std::vector<TString> namesc, std::vector<TString> cuts,
			   TString xname, TString yname, int nbins, float xmin, float xmax, bool uselog, bool norm) {

  std::vector<unsigned int> colors;
  colors.push_back(1);
  colors.push_back(2);
  colors.push_back(4);

  std::vector<unsigned int> lineStyles;
  lineStyles.push_back(1);
  lineStyles.push_back(2);
  lineStyles.push_back(3);

  std::vector<std::vector<TH1F*>> histos; histos.clear();
  for (int itree=0; itree<trees.size(); ++itree) {

    std::vector<TH1F*> histos_tmp; histos_tmp.clear();
    for (int icut=0; icut<cuts.size(); ++icut) {
      TH1F *h_tmp = create1Dhisto(namest[itree],trees[itree],lumi,cuts[icut],var,nbins,xmin,xmax,uselog,colors[itree],lineStyles[itree],"h_"+namest[itree]+"_"+namesc[icut]+"_"+c_name,norm,isdata[itree]);
      h_tmp->SetFillColor(0);
      histos_tmp.push_back(h_tmp);
    }
    histos.push_back(histos_tmp);
  }

  // calculate ratio wrt "0" histogram
  std::cout << " histos = " << histos.size() << "\n";
  std::vector<std::vector<TH1F*>> ratios; ratios.clear();
  for (int ihistos=0; ihistos<histos.size(); ++ihistos) {
    
    std::vector<TH1F*> ratios_tmp; ratios_tmp.clear();
    for (int ihisto=0; ihisto<histos[ihistos].size(); ++ihistos) {
      if (ihisto==0) { continue; }
      TString name = "h_r_"+(TString)histos[ihistos][ihisto]->GetName();
      TH1F* h_r = (TH1F*)histos[ihistos][ihisto]->Clone(name); 
      h_r->Divide(histos[ihistos][ihisto],histos[ihistos][0]);
      h_r->SetLineColor(colors[ihisto-1]); h_r->SetLineStyle(lineStyles[ihisto-1]);
      ratios_tmp.push_back(h_r);
    }
    ratios.push_back(ratios_tmp);
  }
  std::cout << ratios.size() << "\n";
  // plot them
  TCanvas *c_r = new TCanvas("c_ratioOf1Dhistos_"+c_name,"c_ratioOf1Dhistos_"+c_name,500,500);
  for (int ihistos=0; ihistos<ratios.size(); ++ihistos) {
    for (int ihisto=0; ihisto<ratios[ihistos].size(); ++ihisto) {
      
      if ( (ihistos==0) && (ihisto==0) ) {
	ratios[ihistos][ihisto]->GetXaxis()->SetTitle(xname);
	ratios[ihistos][ihisto]->GetYaxis()->SetTitle(yname);
	//ratios[ihistos][ihisto]->GetYaxis()->SetRangeUser(0.,2.);
	ratios[ihistos][ihisto]->Draw("HIST E0");
      }
      else { ratios[ihistos][ihisto]->Draw("HIST E0 sames"); }
    }
  }

  c_r->Print("/uscms_data/d3/loukas/ParticleNetInData/CMSSW_10_2_18/src/PhysicsTools/NanoHRTTools/plotting/plots/particlenet/c_ratioOf1Dhistos_"+c_name+".pdf");
  return c_r;
}





TH1F *create1Dhisto(TString sample_,TTree *tree,TString intLumi,TString cuts,TString branch,int bins,float xmin,float xmax,
                    bool useLog,int color, int style,TString name,bool norm,bool data) {
  TH1::SetDefaultSumw2(kTRUE);

  TString sample;
  if (sample_.Contains("signal"))     { sample = "signal"; }
  if (sample_.Contains("znunu"))      { sample = "qcd"; }
  if (sample_.Contains("znunu-rwgt")) { sample = "qcd-rwgt"; }

  TString cut;
  if (sample == "photons") {
    if (data) { cut ="(passPhoton165_HE10 && "+cuts+")"; }
    else      { cut ="(xsecWeight*puWeight*"+intLumi+")*(passPhoton165_HE10 &&"+cuts+")"; }
  }

  /*  if (sample == "qcd") {
    if (data) { cut ="(passHTTrig && "+cuts+")"; }
    else      { cut ="(xsecWeight*puWeight*"+intLumi+")*(passHTTrig &&"+cuts+")"; }
    }*/


  if (sample == "qcd") {
    if (data) { cut ="("+cuts+")"; }
    else      { cut ="(xsecWeight*puWeight*"+intLumi+")*("+cuts+")"; }                                                                                                             
  }

  if (sample == "qcd-rwgt") {
    if (data) { cut ="("+cuts+")"; }
    else      { cut ="(rewgtfunc(fj_1_tau21)*xsecWeight*puWeight*"+intLumi+")*("+cuts+")"; }
    //else      { cut ="(rewgtfunc(fj_1_sdmass)*xsecWeight*puWeight*"+intLumi+")*("+cuts+")"; }
    //else      { cut ="(rewgtfunc(sqrt((fj_1_sj1_eta-fj_1_sj2_eta)*(fj_1_sj1_eta-fj_1_sj2_eta) + (fj_1_sj1_eta-fj_1_sj2_phi)*(fj_1_sj1_eta-fj_1_sj2_phi)))*xsecWeight*puWeight*"+intLumi+")*("+cuts+")"; }                                                                                                             
  }

  if (sample == "signal") {
    if (data) { cut ="("+cuts+")"; }
    else      { cut ="(xsecWeight*puWeight*"+intLumi+")*("+cuts+"&& fj_1_isH<0.6)"; }
  }

  std::cout << " cut = " << cut << "\n";

  TH1F *hTemp = new TH1F(name,name,bins,xmin,xmax);
  tree->Project(name,branch,cut);

  hTemp->SetLineWidth(3);
  hTemp->SetMarkerSize(0);
  hTemp->SetLineColor(color);
  hTemp->SetFillColor(color);
  hTemp->SetLineStyle(style);

  double error =0.; double integral = hTemp->IntegralAndError(bins,bins+1,error);
  hTemp->SetBinContent(bins,integral);
  hTemp->SetBinError(bins,error);

  if (norm) { hTemp->Scale(1./(hTemp->Integral())); }

  return hTemp;
} //~ end of create1Dhisto




