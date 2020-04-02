/*
 neuron/FalkenuburgerJGP2010.hpp - KCNQ regulation with M1 activation (FalkenburgerJGP2010)

 Copyright (C) 2014 Pranav Kulkarni, Collins Assisi Lab, IISER, Pune <pranavcode@gmail.com>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef INCLUDED_N_FA_HPP
#define INCLUDED_N_FA_HPP

#include "insilico/core/engine.hpp"
#include <random>


namespace insilico {
    
        

class N_F_ER_poisson_spine: public Neuron{
 public:
  void ode_set(state_type &variables, state_type &dxdt, const double t,const unsigned index) {

    double Ek = -87.5, Eampa = 0, Enmda = 0, gamma= 0.06 ; //mV
    //double gp = 0.05; //nS
    double t_p2 = 0.0002, t_p1 = 0.002, tau_s = 0.0043 ; //sec
    double Pw1 = 0.1 , Pw2 = 0.00001 , Pw3 = 3 , Pw4 = 1 ;
    double Mg_spine = 1, eta = 0.33; //mM
    double vm = -20, km=5, vh = -65, kh = -7, tau_m = 0.00008, tau_h = 0.3;
    double Faraday = 96485, R = 8.314, temp = 303.5 ; // C/mol, J/(mol*K), K
    double gammap = 1600, gammad = 300 ;
    double tau_w = 100, P_1 = 1, P_2 = 0.1, P_3 = 0.00001, P_4 = 0 , beta_1 = 80 , beta_2 = 80 , alpha_1 = 0.35 , alpha_2 = 0.55;


    unsigned m_V_index = engine::neuron_index(index,"m_V");
    unsigned h_V_index = engine::neuron_index(index,"h_V");
    unsigned Gbeta_M_index = engine::neuron_index(index,"Gbeta_M");
    unsigned GaGTP_PLC_M_index = engine::neuron_index(index,"GaGTP_PLC_M");
    //unsigned SK_M_index = engine::neuron_index(index,"SK_M");
    unsigned RL_M_index = engine::neuron_index(index,"RL_M");
    unsigned R_M_index = engine::neuron_index(index,"R_M");
    unsigned weight_index = engine::neuron_index(index,"weight");
    unsigned RLG_GDP_M_index = engine::neuron_index(index,"RLG_GDP_M");
    unsigned RLGbeta_M_index = engine::neuron_index(index,"RLGbeta_M");
    unsigned IP3_PH_YFP_C_index = engine::neuron_index(index,"IP3_PH_YFP_C");
    unsigned oxoM_EX_index = engine::neuron_index(index,"oxoM_EX");
    unsigned PH_CFP_PIP2_M_index = engine::neuron_index(index,"PH_CFP_PIP2_M");
    unsigned PH_CFP_C_index = engine::neuron_index(index,"PH_CFP_C");
    unsigned PI4P_M_index = engine::neuron_index(index,"PI4P_M");
    unsigned RGbeta_M_index = engine::neuron_index(index,"RGbeta_M");
    unsigned Ga_GDP_M_index = engine::neuron_index(index,"Ga_GDP_M");
    unsigned PH_YFP_PIP2_M_index = engine::neuron_index(index,"PH_YFP_PIP2_M");
    unsigned DAG_M_index = engine::neuron_index(index,"DAG_M");
    unsigned IP3_C_index = engine::neuron_index(index,"IP3_C");
    unsigned GaGTP_M_index = engine::neuron_index(index,"GaGTP_M");
    unsigned PIP2_M_index = engine::neuron_index(index,"PIP2_M");
    unsigned V_spine_index = engine::neuron_index(index,"V_spine");
    unsigned s_index = engine::neuron_index(index,"s");

    unsigned IP3_D_Cytosol_index = engine::neuron_index(index,"IP3_D_Cytosol");
    unsigned h_ERM_index = engine::neuron_index(index,"h_ERM");

    unsigned h_D_ERM_index = engine::neuron_index(index,"h_D_ERM");
    unsigned Ca_D_Cytosol_index = engine::neuron_index(index,"Ca_D_Cytosol");
    unsigned IP3_Cytosol_index = engine::neuron_index(index,"IP3_Cytosol");
    //unsigned CG_D_Cytosol_index = engine::neuron_index(index,"CG_D_Cytosol");
    unsigned D28kB_high_Cytosol_index = engine::neuron_index(index,"D28kB_high_Cytosol");
    unsigned D28kB_Cytosol_index = engine::neuron_index(index,"D28kB_Cytosol");
    unsigned D28kB_D_Cytosol_index = engine::neuron_index(index,"D28kB_D_Cytosol");
    unsigned D28k_D_Cytosol_index = engine::neuron_index(index,"D28k_D_Cytosol");
    unsigned D28kB_high_D_Cytosol_index = engine::neuron_index(index,"D28kB_high_D_Cytosol");
    unsigned D28k_high_Cytosol_index = engine::neuron_index(index,"D28k_high_Cytosol");
    unsigned D28k_Cytosol_index = engine::neuron_index(index,"D28k_Cytosol");
    unsigned D28k_high_D_Cytosol_index = engine::neuron_index(index,"D28k_high_D_Cytosol");
    //unsigned CG_Cytosol_index = engine::neuron_index(index,"CG_Cytosol");
    //unsigned CGB_D_Cytosol_index = engine::neuron_index(index,"CGB_D_Cytosol");
    unsigned Ca_Cytosol_index = engine::neuron_index(index,"Ca_Cytosol");
    //unsigned PABMg_D_Cytosol_index = engine::neuron_index(index,"PABMg_D_Cytosol");
    //unsigned PA_D_Cytosol_index = engine::neuron_index(index,"PA_D_Cytosol");
    //unsigned PABCa_Cytosol_index = engine::neuron_index(index,"PABCa_Cytosol");
    //unsigned PABMg_Cytosol_index = engine::neuron_index(index,"PABMg_Cytosol");
    //unsigned PABCa_D_Cytosol_index = engine::neuron_index(index,"PABCa_D_Cytosol");
    //unsigned PA_Cytosol_index = engine::neuron_index(index,"PA_Cytosol");
    //unsigned Lpo_index = engine::neuron_index(index,"Lpo");
    //unsigned CGB_Cytosol_index = engine::neuron_index(index,"CGB_Cytosol");


    double s = variables[s_index];
    double V_spine = variables[V_spine_index]; //mV

    double m_V = variables[m_V_index];
    double h_V = variables[h_V_index];    
    double Gbeta_M = variables[Gbeta_M_index];
    double GaGTP_PLC_M = variables[GaGTP_PLC_M_index];
    //double SK_M = variables[SK_M_index];
    double RL_M = variables[RL_M_index];
    double R_M = variables[R_M_index];
    double RLG_GDP_M = variables[RLG_GDP_M_index];
    double RLGbeta_M = variables[RLGbeta_M_index];
    double IP3_PH_YFP_C = variables[IP3_PH_YFP_C_index];
    double oxoM_EX = variables[oxoM_EX_index];
    double PH_CFP_PIP2_M = variables[PH_CFP_PIP2_M_index];
    double PH_CFP_C = variables[PH_CFP_C_index];
    double PI4P_M = variables[PI4P_M_index];
    double RGbeta_M = variables[RGbeta_M_index];
    double Ga_GDP_M = variables[Ga_GDP_M_index];
    double PH_YFP_PIP2_M = variables[PH_YFP_PIP2_M_index];
    double DAG_M = variables[DAG_M_index];
    double IP3_C = variables[IP3_C_index];
    double GaGTP_M = variables[GaGTP_M_index];
    double PIP2_M = variables[PIP2_M_index];
    double weight = variables[weight_index];
    
    double D28kB_high_Cytosol = variables[D28kB_high_Cytosol_index];
    double IP3_D_Cytosol = variables[IP3_D_Cytosol_index];
    double h_ERM = variables[h_ERM_index];
    //double CGB_Cytosol = variables[CGB_Cytosol_index];
    double h_D_ERM = variables[h_D_ERM_index];
    double Ca_D_Cytosol = variables[Ca_D_Cytosol_index];
    double D28kB_Cytosol = variables[D28kB_Cytosol_index];
    //double PABMg_D_Cytosol = variables[PABMg_D_Cytosol_index];
    double IP3_Cytosol = variables[IP3_Cytosol_index];
    //double PA_D_Cytosol = variables[PA_D_Cytosol_index];
    double D28kB_D_Cytosol = variables[D28kB_D_Cytosol_index];
    //double CG_D_Cytosol = variables[CG_D_Cytosol_index];
    double D28kB_high_D_Cytosol = variables[D28kB_high_D_Cytosol_index];
    //double PABCa_Cytosol = variables[PABCa_Cytosol_index];
    double D28k_high_Cytosol = variables[D28k_high_Cytosol_index];
    double D28k_high_D_Cytosol = variables[D28k_high_D_Cytosol_index];
    //double CG_Cytosol = variables[CG_Cytosol_index];
    //double PABMg_Cytosol = variables[PABMg_Cytosol_index];
    //double CGB_D_Cytosol = variables[CGB_D_Cytosol_index];
    double D28k_Cytosol = variables[D28k_Cytosol_index];
    //double PABCa_D_Cytosol = variables[PABCa_D_Cytosol_index];
    double Ca_Cytosol = variables[Ca_Cytosol_index];
    //double PA_Cytosol = variables[PA_Cytosol_index];
    double D28k_D_Cytosol = variables[D28k_D_Cytosol_index];

    // Load Parameters
    double g_Lvdcc = engine::neuron_value(index,"g_Lvdcc");
    double gp = engine::neuron_value(index,"gp");
    double ach_factor = engine::neuron_value(index,"ach_factor");
    double tau1 = engine::neuron_value(index,"tau1") ;
    double tau2 = engine::neuron_value(index,"tau2") ;
    double diff_factor = engine::neuron_value(index,"diff_factor");
    double gn=engine::neuron_value(index,"gn");
    double t_spike = engine::neuron_value(index,"t_spike");
    double g_SK = engine::neuron_value(index,"g_SK"); //nS
    double no_SK = engine::neuron_value(index,"no_SK");
    double g_coupling = engine::neuron_value(index,"g_coupling"); //nS
    double C_spine = engine::neuron_value(index,"C_spine"); //nF
    double V_dendrite = engine::neuron_value(index,"V_dendrite"); //mV
    double carrierValence = engine::neuron_value(index,"carrierValence");
    double Popen = engine::neuron_value(index,"Popen");
    double K4_recover = engine::neuron_value(index,"K4_recover");
    double Kf_reconstitution = engine::neuron_value(index,"Kf_reconstitution");
    double Hill_binding = engine::neuron_value(index,"Hill_binding");
    double K_C_init_uM = engine::neuron_value(index,"K_C_init_uM");
    double netValence_PIP5K_5Pase = engine::neuron_value(index,"netValence_PIP5K_5Pase");
    double KCNQ_PIP2_M_init_molecules_um_2 = engine::neuron_value(index,"KCNQ_PIP2_M_init_molecules_um_2");
    double IP3_PH_YFP_C_init_uM = engine::neuron_value(index,"IP3_PH_YFP_C_init_uM");
    double PH_YFP_PIP2_M_init_molecules_um_2 = engine::neuron_value(index,"PH_YFP_PIP2_M_init_molecules_um_2");
    double KA_PIP2_KCNQ = engine::neuron_value(index,"KA_PIP2_KCNQ");
    double GaGTP_PLC_M_init_molecules_um_2 = engine::neuron_value(index,"GaGTP_PLC_M_init_molecules_um_2");
    double netValence_PLCdiss = engine::neuron_value(index,"netValence_PLCdiss");
    double IP3_C_init_uM = engine::neuron_value(index,"IP3_C_init_uM");
    double Kf_NE_RLG = engine::neuron_value(index,"Kf_NE_RLG");
    double oxoM_EX_init_uM = engine::neuron_value(index,"oxoM_EX_init_uM");
    double GaGTP_M_init_molecules_um_2 = engine::neuron_value(index,"GaGTP_M_init_molecules_um_2");
    double RGbeta_M_init_molecules_um_2 = engine::neuron_value(index,"RGbeta_M_init_molecules_um_2");
    double LumpedJ = engine::neuron_value(index,"LumpedJ");
    double Kr_NE_RG = engine::neuron_value(index,"Kr_NE_RG");
    double netValence_G1beta = engine::neuron_value(index,"netValence_G1beta");
    double Hill_gating = engine::neuron_value(index,"Hill_gating");
    double netValence_PH_YFP_PIP2 = engine::neuron_value(index,"netValence_PH_YFP_PIP2");
    double netValence_L2beta = engine::neuron_value(index,"netValence_L2beta");
    double k5K_rest = engine::neuron_value(index,"k5K_rest");
    double Kf_PLCassoc = engine::neuron_value(index,"Kf_PLCassoc");
    double netValence_PIP2bindKCNQ = engine::neuron_value(index,"netValence_PIP2bindKCNQ");
    double netValence_GTPase_Ga = engine::neuron_value(index,"netValence_GTPase_Ga");
    double speed_PH_IP3 = engine::neuron_value(index,"speed_PH_IP3");
    double Kr_PI4K_4Pase = engine::neuron_value(index,"Kr_PI4K_4Pase");
    double KrG2 = engine::neuron_value(index,"KrG2");
    double netValence_PLC_on_PI4P = engine::neuron_value(index,"netValence_PLC_on_PI4P");
    double Kr_PIP5K_5Pase = engine::neuron_value(index,"Kr_PIP5K_5Pase");
    double alpha = engine::neuron_value(index,"alpha");
    double RG_GDP_M_init_molecules_um_2 = engine::neuron_value(index,"RG_GDP_M_init_molecules_um_2");
    double K_plc = engine::neuron_value(index,"K_plc");
    double GaGDP_PLC_M_init_molecules_um_2 = engine::neuron_value(index,"GaGDP_PLC_M_init_molecules_um_2");
    double DAG_M_init_molecules_um_2 = engine::neuron_value(index,"DAG_M_init_molecules_um_2");
    double IP3_PH_CFP_C_init_uM = engine::neuron_value(index,"IP3_PH_CFP_C_init_uM");
    double Kr_PLCassoc = engine::neuron_value(index,"Kr_PLCassoc");
    double KfG2 = engine::neuron_value(index,"KfG2");
    double PI_M_init_molecules_um_2 = engine::neuron_value(index,"PI_M_init_molecules_um_2");
    double drivingf = engine::neuron_value(index,"drivingf");
    double conductance = engine::neuron_value(index,"conductance");
    double PLC_efficiency_PIP = engine::neuron_value(index,"PLC_efficiency_PIP");
    double k5k_Oxo = engine::neuron_value(index,"k5k_Oxo");
    double Kf_NE_RG = engine::neuron_value(index,"Kf_NE_RG");
    double Size_EX = engine::neuron_value(index,"Size_EX");
    double netValence_L2 = engine::neuron_value(index,"netValence_L2");
    double netValence_L1 = engine::neuron_value(index,"netValence_L1");
    double KL1 = engine::neuron_value(index,"KL1");
    double KCNQ_M_init_molecules_um_2 = engine::neuron_value(index,"KCNQ_M_init_molecules_um_2");
    double Voltage_M = engine::neuron_value(index,"Voltage_M");
    double Ga_GDP_M_init_molecules_um_2 = engine::neuron_value(index,"Ga_GDP_M_init_molecules_um_2");
    double RL_M_init_molecules_um_2 = engine::neuron_value(index,"RL_M_init_molecules_um_2");
    //double mlabfix_F_nmol_ = engine::neuron_value(index,"mlabfix_F_nmol_");
    //double mlabfix_T_ = engine::neuron_value(index,"mlabfix_T_");
    double netValence_NE_RG = engine::neuron_value(index,"netValence_NE_RG");
    //double K_millivolts_per_volt = engine::neuron_value(index,"K_millivolts_per_volt");
    double PLC_basal = engine::neuron_value(index,"PLC_basal");
    double netValence_PLCassoc = engine::neuron_value(index,"netValence_PLCassoc");
    double netValence_NE_GaP = engine::neuron_value(index,"netValence_NE_GaP");
    double RLGbeta_M_init_molecules_um_2 = engine::neuron_value(index,"RLGbeta_M_init_molecules_um_2");
    double netValence_PH_CFP_PIP2 = engine::neuron_value(index,"netValence_PH_CFP_PIP2");
    double K4_rest = engine::neuron_value(index,"K4_rest");
    double PH_YFP_C_init_uM = engine::neuron_value(index,"PH_YFP_C_init_uM");
    double PIP2_M_init_molecules_um_2 = engine::neuron_value(index,"PIP2_M_init_molecules_um_2");
    double PH_CFP_PIP2_M_init_molecules_um_2 = engine::neuron_value(index,"PH_CFP_PIP2_M_init_molecules_um_2");
    double Kr_PLCdiss = engine::neuron_value(index,"Kr_PLCdiss");
    double netValence_PIP2hydr = engine::neuron_value(index,"netValence_PIP2hydr");
    double KD_PH_PIP2 = engine::neuron_value(index,"KD_PH_PIP2");
    double Gbeta_M_init_molecules_um_2 = engine::neuron_value(index,"Gbeta_M_init_molecules_um_2");
    double K_EX_init_uM = engine::neuron_value(index,"K_EX_init_uM");
    double Kf_NE_G = engine::neuron_value(index,"Kf_NE_G");
    //double mlabfix_PI_ = engine::neuron_value(index,"mlabfix_PI_");
    double PI4P_M_init_molecules_um_2 = engine::neuron_value(index,"PI4P_M_init_molecules_um_2");
    //double mlabfix_F_ = engine::neuron_value(index,"mlabfix_F_");
    double netValence_NE_RLG = engine::neuron_value(index,"netValence_NE_RLG");
    double RLG_GDP_M_init_molecules_um_2 = engine::neuron_value(index,"RLG_GDP_M_init_molecules_um_2");
    //double mlabfix_R_ = engine::neuron_value(index,"mlabfix_R_");
    double K_IP3ase = engine::neuron_value(index,"K_IP3ase");
    double Kf_GTPase_Ga = engine::neuron_value(index,"Kf_GTPase_Ga");
    double Kf_PLCdiss = engine::neuron_value(index,"Kf_PLCdiss");
    double G_GDP_M_init_molecules_um_2 = engine::neuron_value(index,"G_GDP_M_init_molecules_um_2");
    double Kr_NE_G = engine::neuron_value(index,"Kr_NE_G");
    //double mlabfix_K_GHK_ = engine::neuron_value(index,"mlabfix_K_GHK_");
    double surface = engine::neuron_value(index,"surface");
    double speed_PIP2_SK = engine::neuron_value(index,"speed_PIP2_SK");
    //double mlabfix_N_pmol_ = engine::neuron_value(index,"mlabfix_N_pmol_");
    double R_M_init_molecules_um_2 = engine::neuron_value(index,"R_M_init_molecules_um_2");
    double KD_PH_IP3 = engine::neuron_value(index,"KD_PH_IP3");
    double Size_M = engine::neuron_value(index,"Size_M"); //um^2
    double Size_C = engine::neuron_value(index,"Size_C"); //um^3
    double netValence_G2 = engine::neuron_value(index,"netValence_G2");
    double KrL1 = engine::neuron_value(index,"KrL1");
    double netValence_G1 = engine::neuron_value(index,"netValence_G1");
    double netValence_PI4K_4Pase = engine::neuron_value(index,"netValence_PI4K_4Pase");
    double netValence_GTPase_GaP = engine::neuron_value(index,"netValence_GTPase_GaP");
    double PH_CFP_C_init_uM = engine::neuron_value(index,"PH_CFP_C_init_uM");
    double netValence_G2beta = engine::neuron_value(index,"netValence_G2beta");
    double speed_PH_PIP2 = engine::neuron_value(index,"speed_PH_PIP2");
    double PLC_M_init_molecules_um_2 = engine::neuron_value(index,"PLC_M_init_molecules_um_2");
    double OxoM_conc_EX_init_uM = engine::neuron_value(index,"OxoM_conc_EX_init_uM");
    double netValence_NE_G = engine::neuron_value(index,"netValence_NE_G");
    double netValence_reconstitution = engine::neuron_value(index,"netValence_reconstitution");
    //double KMOLE = engine::neuron_value(index,"KMOLE");
    double K4_Oxo = engine::neuron_value(index,"K4_Oxo");


    unsigned on1 = 0;
    unsigned on2 = 0;
    unsigned off1 = 0;
    unsigned off2 = 0;

    //Constants

    //double Ks = engine::neuron_value(index,"Ks");
    unsigned number_spikes = engine::neuron_value(index,"number_spikes");
    double LTP_t = engine::neuron_value(index,"LTP_t");
    double LTD_t = engine::neuron_value(index,"LTD_t");
    double P = engine::neuron_value(index,"pump_rate");
    double Ca_threshold = engine::neuron_value(index,"Ca_threshold");
    double tau_Hpo = engine::neuron_value(index,"tau_Hpo");
    double Oxo_on = engine::neuron_value(index,"Oxo_on");
    double Oxo_off= engine::neuron_value(index,"Oxo_off");
    double l_star_IP3deg = engine::neuron_value(index,"l_star_IP3deg");
    double PABMg_F = engine::neuron_value(index,"PABMg_F");
    double r_n_PA_deg = engine::neuron_value(index,"r_n_PA_deg");
    double r_neck_Ca_d = engine::neuron_value(index,"r_neck_Ca_d");
    double mlabfix_F_ = engine::neuron_value(index,"mlabfix_F_");
    double r_spine_PABMg_d = engine::neuron_value(index,"r_spine_PABMg_d");
    double delta = engine::neuron_value(index,"delta");
    double Mg_Cytosol_init_uM = engine::neuron_value(index,"Mg_Cytosol_init_uM");
    double IP3_D_Cytosol_init_uM = engine::neuron_value(index,"IP3_D_Cytosol_init_uM");
    double D_PA_deg = engine::neuron_value(index,"D_PA_deg");
    double lc_IP3deg = engine::neuron_value(index,"lc_IP3deg");
    double h_ERM_init_molecules_um_2 = engine::neuron_value(index,"h_ERM_init_molecules_um_2");
    double CGB_F = engine::neuron_value(index,"CGB_F");
    double l_n_CG_deg = engine::neuron_value(index,"l_n_CG_deg");
    double l_star_PA_deg = engine::neuron_value(index,"l_star_PA_deg");
    double Kf_PA_Ca = engine::neuron_value(index,"Kf_PA_Ca");
    double D_D28k_d = engine::neuron_value(index,"D_D28k_d");
    double r_neck_D28k_high_d = engine::neuron_value(index,"r_neck_D28k_high_d");
    double l_star_PABCa_deg = engine::neuron_value(index,"l_star_PABCa_deg");
    double r_spine_D28kB_high_d = engine::neuron_value(index,"r_spine_D28kB_high_d");
    double Kr_CD28k_high = engine::neuron_value(index,"Kr_CD28k_high");
    double D_D28kB_deg = engine::neuron_value(index,"D_D28kB_deg");
    double r_neck_PABCa_d = engine::neuron_value(index,"r_neck_PABCa_d");
    double Mg_D_Cytosol_init_uM = engine::neuron_value(index,"Mg_D_Cytosol_init_uM");
    double Kr_D28kBDbinding = engine::neuron_value(index,"Kr_D28kBDbinding");
    double r_d_PABMg_deg = engine::neuron_value(index,"r_d_PABMg_deg");
    double lc_PA_deg = engine::neuron_value(index,"lc_PA_deg");
    double PABMg_Cytosol_init_uM = engine::neuron_value(index,"PABMg_Cytosol_init_uM");
    double l_star_D28k_deg = engine::neuron_value(index,"l_star_D28k_deg");
    double D28kB_Cytosol_init_uM = engine::neuron_value(index,"D28kB_Cytosol_init_uM");
    double D_PA_d = engine::neuron_value(index,"D_PA_d");
    double r_D_D28kB_high_deg = engine::neuron_value(index,"r_D_D28kB_high_deg");
    double netValence_flux1 = engine::neuron_value(index,"netValence_flux1");
    double netValence_flux0 = engine::neuron_value(index,"netValence_flux0");
    double l_n_D28k_high_deg = engine::neuron_value(index,"l_n_D28k_high_deg");
    double r_neck_CGB_d = engine::neuron_value(index,"r_neck_CGB_d");
    double Ca_D_ER_init_uM = engine::neuron_value(index,"Ca_D_ER_init_uM");
    double l_CG_d = engine::neuron_value(index,"l_CG_d");
    double D_PABMg_deg = engine::neuron_value(index,"D_PABMg_deg");
    double r_d_CG_deg = engine::neuron_value(index,"r_d_CG_deg");
    double Kf_CGbinding = engine::neuron_value(index,"Kf_CGbinding");
    double l_D28k_d = engine::neuron_value(index,"l_D28k_d");
    double lc_D28kB_high_deg = engine::neuron_value(index,"lc_D28kB_high_deg");
    double l_IP3_d = engine::neuron_value(index,"l_IP3_d");
    double l_n_D28kB_deg = engine::neuron_value(index,"l_n_D28kB_deg");
    double r_neck_PABMg_d = engine::neuron_value(index,"r_neck_PABMg_d");
    double Jch = engine::neuron_value(index,"Jch");
    double r_spine_CG_d = engine::neuron_value(index,"r_spine_CG_d");
    double r_n_D28k_deg = engine::neuron_value(index,"r_n_D28k_deg");
    double Kf_CD28k_high = engine::neuron_value(index,"Kf_CD28k_high");
    double l_star_D28kB_deg = engine::neuron_value(index,"l_star_D28kB_deg");
    double r_neck_PA_d = engine::neuron_value(index,"r_neck_PA_d");
    double Ca_F = engine::neuron_value(index,"Ca_F");
    double PABCa_D_Cytosol_init_uM = engine::neuron_value(index,"PABCa_D_Cytosol_init_uM");
    double Kr_PA_Dbinding = engine::neuron_value(index,"Kr_PA_Dbinding");
    double Ca_Extracellular_init_uM = engine::neuron_value(index,"Ca_Extracellular_init_uM");
    double D_CGB_deg = engine::neuron_value(index,"D_CGB_deg");
    double Size_Extracellular = engine::neuron_value(index,"Size_Extracellular");
    double Size_ERM = engine::neuron_value(index,"Size_ERM");
    double CGB_Cytosol_init_uM = engine::neuron_value(index,"CGB_Cytosol_init_uM");
    double netValence_reaction1 = engine::neuron_value(index,"netValence_reaction1");
    double netValence_reaction0 = engine::neuron_value(index,"netValence_reaction0");
    double r_n_D28kB_high_deg = engine::neuron_value(index,"r_n_D28kB_high_deg");
    double l_star_PABMg_deg = engine::neuron_value(index,"l_star_PABMg_deg");
    double D_D28kB_d = engine::neuron_value(index,"D_D28kB_d");
    double netValence_SERCA_flux = engine::neuron_value(index,"netValence_SERCA_flux");
    double Kr_CG_Dbinding = engine::neuron_value(index,"Kr_CG_Dbinding");
    double D28kB_D_Cytosol_init_uM = engine::neuron_value(index,"D28kB_D_Cytosol_init_uM");
    double Jmax2_IP3R_flux = engine::neuron_value(index,"Jmax2_IP3R_flux");
    double l_n_Ca_deg = engine::neuron_value(index,"l_n_Ca_deg");
    double lc_PABCa_deg = engine::neuron_value(index,"lc_PABCa_deg");
    double lc_D28k_deg = engine::neuron_value(index,"lc_D28k_deg");
    double r_n_PABCa_deg = engine::neuron_value(index,"r_n_PABCa_deg");
    double D28k_Cytosol_init_uM = engine::neuron_value(index,"D28k_Cytosol_init_uM");
    double Kr_PA_Mg = engine::neuron_value(index,"Kr_PA_Mg");
    double K_millivolts_per_volt = engine::neuron_value(index,"K_millivolts_per_volt");
    double mlabfix_K_GHK_ = engine::neuron_value(index,"mlabfix_K_GHK_");
    double Size_PM = engine::neuron_value(index,"Size_PM");
    double netValence_IP3R_fluxD = engine::neuron_value(index,"netValence_IP3R_fluxD");
    double D_D28kB_high_deg = engine::neuron_value(index,"D_D28kB_high_deg");
    double D28kB_F = engine::neuron_value(index,"D28kB_F");
    double Ca_Cytosol_init_uM = engine::neuron_value(index,"Ca_Cytosol_init_uM");
    double r_neck_D28k_d = engine::neuron_value(index,"r_neck_D28k_d");
    double l_CGB_d = engine::neuron_value(index,"l_CGB_d");
    double r_spine_D28kB_d = engine::neuron_value(index,"r_spine_D28kB_d");
    double D_D28k_deg = engine::neuron_value(index,"D_D28k_deg");
    double l_D28k_high_d = engine::neuron_value(index,"l_D28k_high_d");
    double Kdegr_IP3_degr = engine::neuron_value(index,"Kdegr_IP3_degr");
    double lc_D28k_high_deg = engine::neuron_value(index,"lc_D28k_high_deg");
    double D_D28k_high_deg = engine::neuron_value(index,"D_D28k_high_deg");
    double mlabfix_F_nmol_ = engine::neuron_value(index,"mlabfix_F_nmol_");
    double PA_F = engine::neuron_value(index,"PA_F");
    double Voltage_ER = engine::neuron_value(index,"Voltage_ER");
    double Kinh_reaction1 = engine::neuron_value(index,"Kinh_reaction1");
    double Kinh_reaction0 = engine::neuron_value(index,"Kinh_reaction0");
    double lc_D28kB_deg = engine::neuron_value(index,"lc_D28kB_deg");
    double netValence_SERCA_fluxD = engine::neuron_value(index,"netValence_SERCA_fluxD");
    double PABCa_Cytosol_init_uM = engine::neuron_value(index,"PABCa_Cytosol_init_uM");
    double r_n_D28kB_deg = engine::neuron_value(index,"r_n_D28kB_deg");
    double PA_Cytosol_init_uM = engine::neuron_value(index,"PA_Cytosol_init_uM");
    double CGB_D_Cytosol_init_uM = engine::neuron_value(index,"CGB_D_Cytosol_init_uM");
    double D_D28kB_high_d = engine::neuron_value(index,"D_D28kB_high_d");
    double p9 = engine::neuron_value(index,"p9");
    double IP3_Cytosol_init_uM = engine::neuron_value(index,"IP3_Cytosol_init_uM");
    double p8 = engine::neuron_value(index,"p8");
    double p7 = engine::neuron_value(index,"p7");
    double l_n_D28kB_high_deg = engine::neuron_value(index,"l_n_D28kB_high_deg");
    double p6 = engine::neuron_value(index,"p6");
    double p5 = engine::neuron_value(index,"p5");
    double vL_ER_leak_flux = engine::neuron_value(index,"vL_ER_leak_flux");
    double p4 = engine::neuron_value(index,"p4");
    double p3 = engine::neuron_value(index,"p3");
    double p2 = engine::neuron_value(index,"p2");
    double l_n_CGB_deg = engine::neuron_value(index,"l_n_CGB_deg");
    double p1 = engine::neuron_value(index,"p1");
    double p0 = engine::neuron_value(index,"p0");
    double netValence_IP3R_flux = engine::neuron_value(index,"netValence_IP3R_flux");
    double vL_ER_leak_fluxD = engine::neuron_value(index,"vL_ER_leak_fluxD");
    double lc_PABMg_deg = engine::neuron_value(index,"lc_PABMg_deg");
    double Kact_IP3R_flux = engine::neuron_value(index,"Kact_IP3R_flux");
    double IP3_CytosolS = engine::neuron_value(index,"IP3_CytosolS");
    double r_n_PABMg_deg = engine::neuron_value(index,"r_n_PABMg_deg");
    double SVR = engine::neuron_value(index,"SVR");
    double D28kB_high_Cytosol_init_uM = engine::neuron_value(index,"D28kB_high_Cytosol_init_uM");
    double IP3_CytosolD = engine::neuron_value(index,"IP3_CytosolD");
    double JchD = engine::neuron_value(index,"JchD");
    double r_n_CG_deg = engine::neuron_value(index,"r_n_CG_deg");
    double netValence_ER_leak_flux = engine::neuron_value(index,"netValence_ER_leak_flux");
    double r_neck_D28kB_d = engine::neuron_value(index,"r_neck_D28kB_d");
    double mlabfix_T_ = engine::neuron_value(index,"mlabfix_T_");
    double vP_SERCA_fluxD = engine::neuron_value(index,"vP_SERCA_fluxD");
    double CG_Cytosol_init_uM = engine::neuron_value(index,"CG_Cytosol_init_uM");
    double Kdegr_IP3_degr1 = engine::neuron_value(index,"Kdegr_IP3_degr1");
    double D_CG_deg = engine::neuron_value(index,"D_CG_deg");
    double Kf_CG_Dbinding = engine::neuron_value(index,"Kf_CG_Dbinding");
    double l_star_CG_deg = engine::neuron_value(index,"l_star_CG_deg");
    double l_star_D28kB_high_deg = engine::neuron_value(index,"l_star_D28kB_high_deg");
    double r_D_D28k_high_deg = engine::neuron_value(index,"r_D_D28k_high_deg");
    double D_CG_d = engine::neuron_value(index,"D_CG_d");
    double D_IP3_d = engine::neuron_value(index,"D_IP3_d");
    double lc_CG_deg = engine::neuron_value(index,"lc_CG_deg");
    double D28kB_high_F = engine::neuron_value(index,"D28kB_high_F");
    double PABMg_D_Cytosol_init_uM = engine::neuron_value(index,"PABMg_D_Cytosol_init_uM");
    double D28kB_high_D_Cytosol_init_uM = engine::neuron_value(index,"D28kB_high_D_Cytosol_init_uM");
    double kP_SERCA_fluxD = engine::neuron_value(index,"kP_SERCA_fluxD");
    double mlabfix_R_ = engine::neuron_value(index,"mlabfix_R_");
    double Rs = engine::neuron_value(index,"Rs");
    double Kf_PA_MgD = engine::neuron_value(index,"Kf_PA_MgD");
    double r_D_Ca_deg = engine::neuron_value(index,"r_D_Ca_deg");
    double l_Ca_d = engine::neuron_value(index,"l_Ca_d");
    double r_D_D28k_deg = engine::neuron_value(index,"r_D_D28k_deg");
    double r_neck_D28kB_high_d = engine::neuron_value(index,"r_neck_D28kB_high_d");
    double l_n_IP3deg = engine::neuron_value(index,"l_n_IP3deg");
    double r_spine_IP3_d = engine::neuron_value(index,"r_spine_IP3_d");
    double r_spine_Ca_d = engine::neuron_value(index,"r_spine_Ca_d");
    double Ca_ER_init_uM = engine::neuron_value(index,"Ca_ER_init_uM");
    double kP_SERCA_flux = engine::neuron_value(index,"kP_SERCA_flux");
    double tau = engine::neuron_value(index,"tau");
    double h_D_ERM_init_molecules_um_2 = engine::neuron_value(index,"h_D_ERM_init_molecules_um_2");
    double r_neck_CG_d = engine::neuron_value(index,"r_neck_CG_d");
    double l_star_D28k_high_deg = engine::neuron_value(index,"l_star_D28k_high_deg");
    double l_n_PA_deg = engine::neuron_value(index,"l_n_PA_deg");
    double ERDensity_D_ERM_init_molecules_um_2 = engine::neuron_value(index,"ERDensity_D_ERM_init_molecules_um_2");
    double D28k_D_Cytosol_init_uM = engine::neuron_value(index,"D28k_D_Cytosol_init_uM");
    double p11 = engine::neuron_value(index,"p11");
    double D28k_high_D_Cytosol_init_uM = engine::neuron_value(index,"D28k_high_D_Cytosol_init_uM");
    double p10 = engine::neuron_value(index,"p10");
    double l_star_CGB_deg = engine::neuron_value(index,"l_star_CGB_deg");
    double Kf_PA_Mg = engine::neuron_value(index,"Kf_PA_Mg");
    double r_n_Ca_deg = engine::neuron_value(index,"r_n_Ca_deg");
    double CG_D_Cytosol_init_uM = engine::neuron_value(index,"CG_D_Cytosol_init_uM");
    double D_Ca_deg = engine::neuron_value(index,"D_Ca_deg");
    double l_star_Ca_deg = engine::neuron_value(index,"l_star_Ca_deg");
    double Kr_CaD28k_med = engine::neuron_value(index,"Kr_CaD28k_med");
    double Kr_PA_Ca = engine::neuron_value(index,"Kr_PA_Ca");
    double D_CGB_d = engine::neuron_value(index,"D_CGB_d");
    double l_PA_d = engine::neuron_value(index,"l_PA_d");
    double PA_D_Cytosol_init_uM = engine::neuron_value(index,"PA_D_Cytosol_init_uM");
    double l_n_D28k_deg = engine::neuron_value(index,"l_n_D28k_deg");
    double Kf_CD28k_highDbinding = engine::neuron_value(index,"Kf_CD28k_highDbinding");
    double lc_Ca_deg = engine::neuron_value(index,"lc_Ca_deg");
    double Kr_CD28k_highDbinding = engine::neuron_value(index,"Kr_CD28k_highDbinding");
    double Kon_reaction1 = engine::neuron_value(index,"Kon_reaction1");
    double Kon_reaction0 = engine::neuron_value(index,"Kon_reaction0");
    double Size_ER = engine::neuron_value(index,"Size_ER");
    double l_n_PABCa_deg = engine::neuron_value(index,"l_n_PABCa_deg");
    double r_spine_PA_d = engine::neuron_value(index,"r_spine_PA_d");
    double Kf_PA_Dbinding = engine::neuron_value(index,"Kf_PA_Dbinding");
    double D28k_high_F = engine::neuron_value(index,"D28k_high_F");
    double Ca_D_Cytosol_init_uM = engine::neuron_value(index,"Ca_D_Cytosol_init_uM");
    double r_n_CGB_deg = engine::neuron_value(index,"r_n_CGB_deg");
    double r_spine_CGB_d = engine::neuron_value(index,"r_spine_CGB_d");
    double CG_F = engine::neuron_value(index,"CG_F");
    double r_n_D28k_high_deg = engine::neuron_value(index,"r_n_D28k_high_deg");
    double ERDensity_ERM_init_molecules_um_2 = engine::neuron_value(index,"ERDensity_ERM_init_molecules_um_2");
    double vP_SERCA_flux = engine::neuron_value(index,"vP_SERCA_flux");
    double D28k_F = engine::neuron_value(index,"D28k_F");
    double Kf_CaD28k_med = engine::neuron_value(index,"Kf_CaD28k_med");
    double D28k_high_Cytosol_init_uM = engine::neuron_value(index,"D28k_high_Cytosol_init_uM");
    double mlabfix_PI_ = engine::neuron_value(index,"mlabfix_PI_");
    double r_spine_D28k_d = engine::neuron_value(index,"r_spine_D28k_d");
    double l_PABCa_d = engine::neuron_value(index,"l_PABCa_d");
    double Kr_CGbinding = engine::neuron_value(index,"Kr_CGbinding");
    double Ks = engine::neuron_value(index,"Ks");
    double l_D28kB_high_d = engine::neuron_value(index,"l_D28kB_high_d");
    double IP3_F = engine::neuron_value(index,"IP3_F");
    double D_PABCa_d = engine::neuron_value(index,"D_PABCa_d");
    double Js = engine::neuron_value(index,"Js");
    double Size_Cytosol = engine::neuron_value(index,"Size_Cytosol");
    double r_spine_D28k_high_d = engine::neuron_value(index,"r_spine_D28k_high_d");
    double mlabfix_N_pmol_ = engine::neuron_value(index,"mlabfix_N_pmol_");
    double Kf_D28kBDbinding = engine::neuron_value(index,"Kf_D28kBDbinding");
    double l_PABMg_d = engine::neuron_value(index,"l_PABMg_d");
    double l_n_PABMg_deg = engine::neuron_value(index,"l_n_PABMg_deg");
    double l_D28kB_d = engine::neuron_value(index,"l_D28kB_d");
    double t2_flux1 = engine::neuron_value(index,"t2_flux1");
    double t2_flux0 = engine::neuron_value(index,"t2_flux0");
    double r_d_PABCa_deg = engine::neuron_value(index,"r_d_PABCa_deg");
    double Kact_IP3R_fluxD = engine::neuron_value(index,"Kact_IP3R_fluxD");
    double r_neck_IP3_d = engine::neuron_value(index,"r_neck_IP3_d");
    double KMOLE = engine::neuron_value(index,"KMOLE");
    double D_PABMg_d = engine::neuron_value(index,"D_PABMg_d");
    double Voltage_PM = engine::neuron_value(index,"Voltage_PM");
    double dI_IP3R_flux = engine::neuron_value(index,"dI_IP3R_flux");
    double D_Ca_d = engine::neuron_value(index,"D_Ca_d");
    double D_D28k_high_d = engine::neuron_value(index,"D_D28k_high_d");
    double PABCa_F = engine::neuron_value(index,"PABCa_F");
    double Ca_D_Extracellular_init_uM = engine::neuron_value(index,"Ca_D_Extracellular_init_uM");
    double r_d_CGB_deg = engine::neuron_value(index,"r_d_CGB_deg");
    double r_D_D28kB_deg = engine::neuron_value(index,"r_D_D28kB_deg");
    double r_spine_PABCa_d = engine::neuron_value(index,"r_spine_PABCa_d");
    double r_d_IP3deg = engine::neuron_value(index,"r_d_IP3deg");
    double D_PABCa_deg = engine::neuron_value(index,"D_PABCa_deg");
    double r_n_IP3deg = engine::neuron_value(index,"r_n_IP3deg");
    double Jmax2_IP3R_fluxD = engine::neuron_value(index,"Jmax2_IP3R_fluxD");
    double Kr_PA_MgD = engine::neuron_value(index,"Kr_PA_MgD");
    double netValence_ER_leak_fluxD = engine::neuron_value(index,"netValence_ER_leak_fluxD");
    double r_d_PA_deg = engine::neuron_value(index,"r_d_PA_deg");
    double D_IP3deg = engine::neuron_value(index,"D_IP3deg");
    double t1_flux1 = engine::neuron_value(index,"t1_flux1");
    double dI_IP3R_fluxD = engine::neuron_value(index,"dI_IP3R_fluxD");
    double t1_flux0 = engine::neuron_value(index,"t1_flux0");
    double lc_CGB_deg = engine::neuron_value(index,"lc_CGB_deg");
    double frequency = engine::neuron_value(index,"frequency");
    double J_IP3=0;
/*    for(int idx=0; idx<=11; idx++){

        double p = engine::neuron_value(index,"p"+std::to_string(idx));
        short flg = 0;
        if (t > idx*tau + delta){
            flg = 1;
        }
        else{
            flg = 0;
        }
        double J_pulse = p*flg*exp(-1*((t - (idx*tau +delta))*Ks)) ;
        
        J_IP3 = J_IP3 + J_pulse;
    }
*/

    if (t>Oxo_on){
      on1 = 1;
      off1 = 0;
    }
    else{
      on1 = 0;
      off1 = 1;
    }
    if (t<Oxo_off){
      on2 = 1;
      off2 = 0;
    }
    else{
      on2 = 0;
      off2 = 1;
    }
    if(on1 && on2){
        
        engine::neuron_value(index,"V_dendrite",-71.50);
    }
    else{
        
        engine::neuron_value(index,"V_dendrite",-75.00);        
    }
    // Functions for calculations
    double period = 1/frequency;
    double Ca_eq =  (on1 && on2) ? 0.150 : Ca_D_Cytosol; // uM
    double KCNQ_PIP2_M = 0;
    double Eca =  (on1 && on2) ? 1.420 : 0.320; // uM
    double EC_PIP2 = (on1 && on2) ? 3010 : 365; // molecules/um2
    double UnitFactor_uM_um3_molecules_neg_1 = (1.0 / 602.0);
    double K_OxoM_conc_EX_total = ((Size_EX * OxoM_conc_EX_init_uM) + (Size_EX * oxoM_EX_init_uM) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLG_GDP_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RL_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLGbeta_M_init_molecules_um_2));
    double Kf_PH_YFP_IP3 = speed_PH_IP3;
    double Kr_PH_YFP_IP3 = (KD_PH_IP3 * speed_PH_IP3);
    double K_PH_YFP_C_total = ((Size_C * PH_YFP_C_init_uM) + (Size_C * IP3_PH_YFP_C_init_uM) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * PH_YFP_PIP2_M_init_molecules_um_2));
    double PH_YFP_C = ((K_PH_YFP_C_total - (Size_C * IP3_PH_YFP_C) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * PH_YFP_PIP2_M)) / Size_C);
    double J_PH_YFP_IP3 = (((Kf_PH_YFP_IP3 * IP3_C) * PH_YFP_C) - (Kr_PH_YFP_IP3 * IP3_PH_YFP_C));
    double KfL1 = (KrL1 / KL1);
    double KfL2 = KfL1;
    double K_EX = K_EX_init_uM;
    double Kr_PH_YFP_PIP2 = (KD_PH_PIP2 * speed_PH_PIP2);
    double Kf_PH_YFP_PIP2 = speed_PH_PIP2;
    double J_PH_YFP_PIP2 = (((Kf_PH_YFP_PIP2 * PH_YFP_C) * PIP2_M) - (Kr_PH_YFP_PIP2 * PH_YFP_PIP2_M));
    double K_KCNQ_PIP2_M_total = ((Size_M * UnitFactor_uM_um3_molecules_neg_1 * KCNQ_PIP2_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * KCNQ_M_init_molecules_um_2));
    //double SK_PIP2_M = (no_SK - SK_M);
    double K_GaGDP_PLC_M_total = ((Size_M * UnitFactor_uM_um3_molecules_neg_1 * GaGDP_PLC_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * Ga_GDP_M_init_molecules_um_2) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RGbeta_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * GaGTP_M_init_molecules_um_2) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLGbeta_M_init_molecules_um_2) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * Gbeta_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * GaGTP_PLC_M_init_molecules_um_2));
    double GaGDP_PLC_M = ((K_GaGDP_PLC_M_total - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * Ga_GDP_M) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RGbeta_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * GaGTP_M) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLGbeta_M) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * Gbeta_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * GaGTP_PLC_M)) / (Size_M * UnitFactor_uM_um3_molecules_neg_1));
    double J_NE_GaP = (4.7 * GaGDP_PLC_M);
    double K_IP3_PH_CFP_C_total = ((Size_C * IP3_PH_CFP_C_init_uM) + (Size_C * PH_CFP_C_init_uM) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * PH_CFP_PIP2_M_init_molecules_um_2));
    double IP3_PH_CFP_C = ((K_IP3_PH_CFP_C_total - (Size_C * PH_CFP_C) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * PH_CFP_PIP2_M)) / Size_C);
    double allPH_cytoplasm = (PH_CFP_C + PH_YFP_C + IP3_PH_CFP_C + IP3_PH_YFP_C);
    double K_G_GDP_M_total = ((Size_M * UnitFactor_uM_um3_molecules_neg_1 * G_GDP_M_init_molecules_um_2) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * R_M_init_molecules_um_2) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RL_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * Gbeta_M_init_molecules_um_2));
    double G_GDP_M = ((K_G_GDP_M_total + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * R_M) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RL_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * Gbeta_M)) / (Size_M * UnitFactor_uM_um3_molecules_neg_1));
    double J_NE_G = ((Kf_NE_G * G_GDP_M) - ((Kr_NE_G * GaGTP_M) * Gbeta_M));
    double J_PLC_on_PI4P = (PI4P_M * PLC_efficiency_PIP * (PLC_basal + (GaGTP_PLC_M * K_plc)));
    double Kr_PH_CFP_IP3 = (KD_PH_IP3 * speed_PH_IP3);
    double Kf_PH_CFP_PIP2 = speed_PH_PIP2;
    double K_PLC_M_total = ((Size_M * UnitFactor_uM_um3_molecules_neg_1 * PLC_M_init_molecules_um_2) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * Ga_GDP_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RGbeta_M_init_molecules_um_2) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * GaGTP_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLGbeta_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * Gbeta_M_init_molecules_um_2));
    double PLC_M = ((K_PLC_M_total + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * Ga_GDP_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RGbeta_M) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * GaGTP_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLGbeta_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * Gbeta_M)) / (Size_M * UnitFactor_uM_um3_molecules_neg_1));
    double J_PLCdiss = ((Kf_PLCdiss * GaGDP_PLC_M) - ((Kr_PLCdiss * PLC_M) * Ga_GDP_M));
    //double Kr_PIP2bindSK = (on1 && on2) ? 0.3 : 0.1;
    //double KD_PIP2_SK = pow(EC_PIP2, Hill_binding);
    //double Kf_PIP2bindSK = (Kr_PIP2bindSK / KD_PIP2_SK);
    //double J_PIP2bindSK = ((Kf_PIP2bindSK * pow(PIP2_M , Hill_binding) * SK_M) - (Kr_PIP2bindSK * SK_PIP2_M ));
    double allRL = (RL_M + RLG_GDP_M + RLGbeta_M);
    double KG2 = (KrG2 / KfG2);
    double KG1 = (KG2 * alpha);
    double KfG1 = (0.1 * KfG2);
    double KrG1 = (KfG1 * KG1);
    double Kr_G1beta = KrG1;
    double J_NE_RLG = (Kf_NE_RLG * RLG_GDP_M);
    double KL2 = (KL1 / alpha);
    double KrL2 = (KL2 * KfL2);
    double Kr_L2beta = KrL2;
    double PI_M = PI_M_init_molecules_um_2;
    double Kf_L2 = KfL2;
    double Kf_L1 = KfL1;
    double Kr_G2 = KrG2;
    double Kr_G1 = KrG1;
    double Kf_Oxo_appl = (50.0 * on1 * on2);
    double K_RG_GDP_M_total = ((Size_M * UnitFactor_uM_um3_molecules_neg_1 * RG_GDP_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLG_GDP_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RGbeta_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * R_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RL_M_init_molecules_um_2) + (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLGbeta_M_init_molecules_um_2));
    double RG_GDP_M = ((K_RG_GDP_M_total - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLG_GDP_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RGbeta_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * R_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RL_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLGbeta_M)) / (Size_M * UnitFactor_uM_um3_molecules_neg_1));
    double Kr_L2 = KrL2;
    double J_L2 = (((Kf_L2 * oxoM_EX) * RG_GDP_M) - (Kr_L2 * RLG_GDP_M));
    double Kr_L1 = KrL1;
    double J_L1 = (((Kf_L1 * R_M) * oxoM_EX) - (Kr_L1 * RL_M));
    double Kf_PH_CFP_IP3 = speed_PH_IP3;
    double Kf_G2beta = KfG2;
    double Kr_PH_CFP_PIP2 = (KD_PH_PIP2 * speed_PH_PIP2);
    double J_PH_CFP_PIP2 = (((Kf_PH_CFP_PIP2 * PIP2_M) * PH_CFP_C) - (Kr_PH_CFP_PIP2 * PH_CFP_PIP2_M));
    double Kr_G2beta = KrG2;
    double J_G2beta = (((Kf_G2beta * RL_M) * Gbeta_M) - (Kr_G2beta * RLGbeta_M));
    double Kf_PI4K_4Pase = ((K4_recover * off2) + (K4_rest * off1) + (K4_Oxo * on1 * on2));
    double J_reconstitution = (Kf_reconstitution * Gbeta_M * Ga_GDP_M);
    double J_PH_CFP_IP3 = (((Kf_PH_CFP_IP3 * IP3_C) * PH_CFP_C) - (Kr_PH_CFP_IP3 * IP3_PH_CFP_C));
    double J_GTPase_GaP = (15.0 * GaGTP_PLC_M);
    double allGaPLC = (GaGTP_PLC_M + GaGDP_PLC_M);
    double Kr_Oxo_appl = (50.0 * off2);
    double J_PI4K_4Pase = ((Kf_PI4K_4Pase * PI_M) - (Kr_PI4K_4Pase * PI4P_M));
    double J_NE_RG = ((Kf_NE_RG * RG_GDP_M) - ((Kr_NE_RG * GaGTP_M) * RGbeta_M));
    double J_IP3Pase = (K_IP3ase * IP3_C);
    double Kf_PIP5K_5Pase = ((off1 * k5K_rest) + (on1 * on2 * k5k_Oxo) + (off2 * k5K_rest));
    double Kf_G1beta = KfG1;
    double Kf_L2beta = KfL2;
    double J_PLCassoc = (((Kf_PLCassoc * PLC_M) * GaGTP_M) - (Kr_PLCassoc * GaGTP_PLC_M));
    double Kf_G2 = KfG2;
    double Kf_G1 = KfG1;
    double KFlux_M_EX = (Size_M / Size_EX);
    double J_G1beta = (((Kf_G1beta * Gbeta_M) * R_M) - (Kr_G1beta * RGbeta_M));
    double J_L2beta = (((Kf_L2beta * RGbeta_M) * oxoM_EX) - (Kr_L2beta * RLGbeta_M));
    double J_PIP2hydr = (PIP2_M * (PLC_basal + (K_plc * GaGTP_PLC_M)));
    double J_G2 = (((Kf_G2 * G_GDP_M) * RL_M) - (Kr_G2 * RLG_GDP_M));
    double J_G1 = (((Kf_G1 * G_GDP_M) * R_M) - (Kr_G1 * RG_GDP_M));
    double J_GTPase_Ga = (Kf_GTPase_Ga * GaGTP_M);
    double OxoM_conc_EX = ((K_OxoM_conc_EX_total - (Size_EX * oxoM_EX) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLG_GDP_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RL_M) - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * RLGbeta_M)) / Size_EX);
    double J_Oxo_appl = ((Kf_Oxo_appl * OxoM_conc_EX) - (Kr_Oxo_appl * oxoM_EX));
    double allPIP2M = (PIP2_M + KCNQ_PIP2_M + PH_CFP_PIP2_M + PH_YFP_PIP2_M);
    double IM = (surface * Popen * conductance * drivingf * pow(KCNQ_PIP2_M , Hill_gating));
    double KFlux_M_C = (Size_M / Size_C);
    double J_PIP5K_5Pase = ((Kf_PIP5K_5Pase * PI4P_M) - (Kr_PIP5K_5Pase * PIP2_M));
    double K_C = K_C_init_uM;
    double allRGbeta = (RG_GDP_M + RLG_GDP_M + RLGbeta_M + RGbeta_M);
    double allGabg = (G_GDP_M + RG_GDP_M + RLG_GDP_M);
    double s_sat = pow(Ca_Cytosol,4.1)/(pow(Ca_Cytosol,4.1) + pow(Eca,4.1)); //Ca in uM
    double J_pulses = SVR*Js*J_IP3;
    double ERDensity_ERM = ERDensity_ERM_init_molecules_um_2;
    double J_SERCA_flux = (ERDensity_ERM * vP_SERCA_flux * Ca_Cytosol * Ca_Cytosol / ((kP_SERCA_flux * kP_SERCA_flux) + (Ca_Cytosol * Ca_Cytosol)));
    //double J_CGB_deg = ((D_CGB_deg * pow(r_n_CGB_deg,2.0) * (CGB_D_Cytosol - CGB_Cytosol) / l_n_CGB_deg / pow(r_d_CGB_deg,2.0) / l_star_CGB_deg) + (D_CGB_deg * (CGB_D_Cytosol - CGB_F) / l_star_CGB_deg / lc_CGB_deg));
    double KFlux_PM_Cytosol = (Size_PM / Size_Cytosol);
    double J_D28kB_high_d = (0.75 * D_D28kB_high_d * (D28kB_high_Cytosol - D28kB_high_D_Cytosol) * pow(r_neck_D28kB_high_d,2.0) / l_D28kB_high_d / pow(r_spine_D28kB_high_d,3.0));
    double J_D28kB_d = (0.75 * D_D28kB_d * (D28kB_Cytosol - D28kB_Cytosol) * pow(r_neck_D28kB_d,2.0) / l_D28kB_d / pow(r_spine_D28kB_d,3.0));
    double J_reaction1 =  - ((Kinh_reaction1 - ((Ca_D_Cytosol + Kinh_reaction1) * h_D_ERM)) * Kon_reaction1);
    double J_reaction0 =  - ((Kinh_reaction0 - ((Ca_Cytosol + Kinh_reaction0) * h_ERM)) * Kon_reaction0);
    double ERDensity_D_ERM = ERDensity_D_ERM_init_molecules_um_2;
    double Ca_D_ER = Ca_D_ER_init_uM;
    double J_IP3R_fluxD =  - (ERDensity_D_ERM * Jmax2_IP3R_fluxD * (1.0 - (Ca_D_Cytosol / Ca_D_ER)) * pow((h_D_ERM * IP3_D_Cytosol * Ca_D_Cytosol / (IP3_D_Cytosol + dI_IP3R_fluxD) / (Ca_D_Cytosol + Kact_IP3R_fluxD)),3.0));
    //double Mg_Cytosol = Mg_Cytosol_init_uM;
    double J_IP3degrade_IP3_degr = (Kdegr_IP3_degr * (IP3_Cytosol - IP3_CytosolS));
    double KFlux_ERM_Cytosol = (Size_ERM / Size_Cytosol);
    double J_D28k_high_d = (0.75 * D_D28k_high_d * (D28k_high_Cytosol - D28k_high_D_Cytosol) * pow(r_neck_D28k_high_d,2.0) / l_D28k_high_d / pow(r_spine_D28k_high_d,3.0));
    double Ca_ER = Ca_ER_init_uM;
    double J_IP3R_flux =  - (ERDensity_ERM * Jmax2_IP3R_flux * (1.0 - (Ca_Cytosol / Ca_ER)) * (pow(h_ERM * IP3_C * Ca_Cytosol / (IP3_C + dI_IP3R_flux) / (Ca_Cytosol + Kact_IP3R_flux),3.0)));
    //double J_PABCa_deg = ((D_PABCa_deg * pow(r_n_PABCa_deg,2.0) * (PABCa_D_Cytosol - PABCa_Cytosol) / l_n_PABCa_deg / pow(r_d_PABCa_deg,2.0) / l_star_PABCa_deg) + (D_PABCa_deg * (PABCa_D_Cytosol - PABCa_F) / l_star_PABCa_deg / lc_PABCa_deg));
    //double J_CG_deg = ((D_CG_deg * pow(r_n_CG_deg,2.0) * (CG_D_Cytosol - CG_Cytosol) / l_n_CG_deg / pow(r_d_CG_deg,2.0) / l_star_CG_deg) + (D_CG_deg * (CG_D_Cytosol - CG_F) / l_star_CG_deg / lc_CG_deg));
    double J_Ca_d = diff_factor*(0.75 * D_Ca_d * (Ca_Cytosol - Ca_eq/*Ca_D_Cytosol*/) * pow(r_neck_Ca_d,2.0) / l_Ca_d / pow(r_spine_Ca_d,3.0));
    double J_D28kB_deg = ((D_D28kB_deg * pow(r_n_D28kB_deg,2.0) * (D28kB_D_Cytosol - D28kB_Cytosol) / l_n_D28kB_deg / pow(r_D_D28kB_deg,2.0) / l_star_D28kB_deg) + (D_D28kB_deg * (D28kB_D_Cytosol - D28kB_F) / l_star_D28kB_deg / lc_D28kB_deg));
    double J_D28k_deg = ((D_D28k_deg * pow(r_n_D28k_deg,2.0) * (D28k_D_Cytosol - D28k_Cytosol) / l_n_D28k_deg / pow(r_D_D28k_deg,2.0) / l_star_D28k_deg) + (D_D28k_deg * (D28k_D_Cytosol - D28k_F) / l_star_D28k_deg / lc_D28k_deg));
    double J_ER_leak_fluxD =  - (ERDensity_D_ERM * vL_ER_leak_fluxD * (1.0 - (Ca_D_Cytosol / Ca_D_ER)));
    double J_D28k_high_deg = ((D_D28k_high_deg * pow(r_n_D28k_high_deg,2.0) * (D28k_high_D_Cytosol - D28k_high_Cytosol) / l_n_D28k_high_deg / pow(r_D_D28k_high_deg,2.0) / l_star_D28k_high_deg) + (D_D28k_high_deg * (D28k_high_D_Cytosol - D28k_high_F) / l_star_D28k_high_deg / lc_D28k_high_deg));
    double J_D28kB_high_deg = ((D_D28kB_high_deg * pow(r_n_D28kB_high_deg,2.0) * (D28kB_high_D_Cytosol - D28kB_high_Cytosol) / l_n_D28kB_high_deg / pow(r_D_D28kB_high_deg,2.0) / l_star_D28kB_high_deg) + (D_D28kB_high_deg * (D28kB_high_D_Cytosol - D28kB_high_F) / l_star_D28kB_high_deg / lc_D28kB_high_deg));
    //double J_PABMg_deg = ((D_PABMg_deg * pow(r_n_PABMg_deg,2.0) * (PABMg_D_Cytosol - PABMg_Cytosol) / l_n_PABMg_deg / pow(r_d_PABMg_deg,2.0) / l_star_PABMg_deg) + (D_PABMg_deg * (PABMg_D_Cytosol - PABMg_F) / l_star_PABMg_deg / lc_PABMg_deg));
    double J_D28kBDbinding = (((Kf_D28kBDbinding * D28k_D_Cytosol) * Ca_D_Cytosol) - (Kr_D28kBDbinding * D28kB_D_Cytosol));
    double Ca_D_Extracellular = Ca_D_Extracellular_init_uM;
    double J_flux1 = 0;//(JchD * (t > t1_flux1) * (t < t2_flux1) * (Ca_D_Extracellular - Ca_D_Cytosol));
    double Ca_Extracellular = Ca_Extracellular_init_uM;
    double J_flux0 = 0; //(Jch * (t > t1_flux0) * (t < t2_flux0) * (Ca_Extracellular - Ca_Cytosol));
    //double J_PA_d = (0.75 * D_PA_d * (PA_Cytosol - PA_D_Cytosol) * pow(r_neck_PA_d,2.0) / l_PA_d / pow(r_spine_PA_d,3.0));
    //double J_PA_Dbinding = (((Kf_PA_Dbinding * PA_D_Cytosol) * Ca_D_Cytosol) - (Kr_PA_Dbinding * PABCa_D_Cytosol));
    double J_Ca_deg = ((D_Ca_deg * pow(r_n_Ca_deg,2.0) * (Ca_D_Cytosol - Ca_Cytosol) / l_n_Ca_deg / pow(r_D_Ca_deg,2.0) / l_star_Ca_deg) + (D_Ca_deg * (Ca_D_Cytosol - Ca_F) / l_star_Ca_deg / lc_Ca_deg));
    //double J_CG_Dbinding = (((Kf_CG_Dbinding * Ca_D_Cytosol) * CG_D_Cytosol) - (Kr_CG_Dbinding * CGB_D_Cytosol));
    double J_CD28k_highDbinding = (((Kf_CD28k_highDbinding * Ca_D_Cytosol) * D28k_high_D_Cytosol) - (Kr_CD28k_highDbinding * D28kB_high_D_Cytosol));
    double J_CaD28k_med = (((Kf_CaD28k_med * D28k_Cytosol) * Ca_Cytosol) - (Kr_CaD28k_med * D28kB_Cytosol));
    //double Mg_D_Cytosol = Mg_D_Cytosol_init_uM;
    //double J_PACa = (((Kf_PA_Ca * PA_Cytosol) * Ca_Cytosol) - (Kr_PA_Ca * PABCa_Cytosol));
    //double J_PABCa_d = (0.75 * D_PABCa_d * (PABCa_Cytosol - PABCa_D_Cytosol) * pow(r_neck_PABCa_d,2.0) / l_PABCa_d / pow(r_spine_PABCa_d,3.0));
    double J_IP3_d = (0.75 * D_IP3_d * (IP3_Cytosol - IP3_D_Cytosol) * pow(r_neck_IP3_d,2.0) / l_IP3_d / pow(r_spine_IP3_d,3.0));
    double J_SERCA_fluxD = (ERDensity_D_ERM * vP_SERCA_fluxD * Ca_D_Cytosol * Ca_D_Cytosol / ((kP_SERCA_fluxD * kP_SERCA_fluxD) + (Ca_D_Cytosol * Ca_D_Cytosol)));
    double J_ER_leak_flux =  - (ERDensity_ERM * vL_ER_leak_flux * (1.0 - (Ca_Cytosol / Ca_ER)));
    //double J_PABMg_d = (0.75 * D_PABMg_d * (PABMg_Cytosol - PABMg_D_Cytosol) * pow(r_neck_PABMg_d,2.0) / l_PABMg_d / pow(r_spine_PABMg_d,3.0));
    double J_IP3degrade_IP3_degr1 = (Kdegr_IP3_degr1 * (IP3_D_Cytosol - IP3_CytosolD));
    //double J_PAMg_PA_MgD = (((Kf_PA_MgD * Mg_D_Cytosol) * PA_D_Cytosol) - (Kr_PA_MgD * PABMg_D_Cytosol));
    double J_IP3deg = ((D_IP3deg * pow(r_n_IP3deg,2.0) * (IP3_D_Cytosol - IP3_Cytosol) / l_n_IP3deg / pow(r_d_IP3deg,2.0) / l_star_IP3deg) + (D_IP3deg * (IP3_D_Cytosol - IP3_F) / l_star_IP3deg / lc_IP3deg));
    //double J_CGB_d = (0.75 * D_CGB_d * (CGB_Cytosol - CGB_D_Cytosol) * pow(r_neck_CGB_d,2.0) / l_CGB_d / pow(r_spine_CGB_d,3.0));
    //double J_PAMg_PA_Mg = (((Kf_PA_Mg * PA_Cytosol) * Mg_Cytosol) - (Kr_PA_Mg * PABMg_Cytosol));
    //double J_PA_deg = ((D_PA_deg * pow(r_n_PA_deg,2.0) * (PA_D_Cytosol - PA_Cytosol) / l_n_PA_deg / pow(r_d_PA_deg,2.0) / l_star_PA_deg) + (D_PA_deg * (PA_D_Cytosol - PA_F) / l_star_PA_deg / lc_PA_deg));
    double J_D28k_d = (0.75 * D_D28k_d * (D28k_Cytosol - D28k_D_Cytosol) * pow(r_neck_D28k_d,2.0) / l_D28k_d / pow(r_spine_D28k_d,3.0));
    //double J_CGbinding = (((Kf_CGbinding * Ca_Cytosol) * CG_Cytosol) - (Kr_CGbinding * CGB_Cytosol));
    //double J_CG_d = (0.75 * D_CG_d * (CG_Cytosol - CG_D_Cytosol) * pow(r_neck_CG_d,2.0) / l_CG_d / pow(r_spine_CG_d,3.0));
    double J_CD28k_high = (((Kf_CD28k_high * Ca_Cytosol) * D28k_high_Cytosol) - (Kr_CD28k_high * D28kB_high_Cytosol));
    double H_SK = pow(Ca_Cytosol,4.1)/(pow(Ca_Cytosol,4.1) + pow(Eca,4.1));
    double H_PIP2 = (pow(PIP2_M,Hill_gating)/(pow(PIP2_M,Hill_gating) + pow(EC_PIP2 ,Hill_gating)));
    double I_dendrite = g_coupling*(V_spine - V_dendrite);

    double m_inf = 1/(1 + exp((vm - V_spine)/km));
    double h_inf = 1/(1 + exp((vh - V_spine)/kh));

    double Theta_Ca = 0.25 + exp(beta_2*(Ca_Cytosol - alpha_2))/(1 + exp(beta_2*(Ca_Cytosol - alpha_2))) + -0.25*exp(beta_1*(Ca_Cytosol - alpha_1))/(1 + exp(beta_1*(Ca_Cytosol - alpha_1))) ;
    P_4 = 3*Theta_Ca ;
    tau_w = Pw4 + Pw1/(Pw2 + pow(Ca_Cytosol,Pw3)) ;
    double I_NMDA = 0 ;
    double I_AMPA = 0 ;
    double gnmda = 0 ;
    double Flux_NMDA = 0 ;
    double plasticity_rate = 0;
    double ltp_on = (Ca_Cytosol > LTP_t) ? 1.0 : 0.0;
    double ltd_on = (Ca_Cytosol > LTD_t) ? 1.0 : 0.0;
    double Pf_NMDA = 0.135 ;
    double ECa = 2.303*8.314*309150*log10(1000/Ca_Cytosol)/(2*96485) ;
    double rn = 0;

    if (t > t_spike){
        
        unsigned counter = 0;

        double t_spikes = t_spike;

        //rn = ((exp((t_spike - t)/tau1) - exp((t_spike - t)/tau2))/(1 + eta*Mg_spine*exp(-1*gamma*V_spine))) ;
        //I_NMDA = gn*rn*(V_spine - Enmda) ;
        //I_AMPA = gp*(exp((t_spike - t)/t_p1) - exp((t_spike - t)/t_p2))*(V_spine - Eampa);
        //Flux_NMDA = (5182)*(Pf_NMDA/4)*rn*gn*(V_spine - ECa) ; //I_NMDA*0.1*25900;
        //t_spikes = t_spikes + 0.01;
        for (auto& it: insilico::engine::periods) { //iterate over spikes


            t_spikes = t_spikes + it ;

            if(t > t_spikes){
                rn = ((exp((t_spikes - t)/tau1) - exp((t_spikes - t)/tau2))/(1 + eta*Mg_spine*exp(-1*gamma*V_spine))) ;
                I_NMDA += gn*rn*(V_spine - Enmda) ;
                I_AMPA += gp*(exp((t_spikes - t)/t_p1) - exp((t_spikes - t)/t_p2))*(V_spine - Eampa);
                Flux_NMDA += (5182)*(Pf_NMDA/4)*rn*gn*(V_spine - ECa) ; //I_NMDA*0.1*25900;

            }
            if (counter< number_spikes){
                 counter +=1 ;
            }



        }
        
    }
    double qca_Nav = 1.927346;
    double I_SK = s*H_PIP2*g_SK*(V_spine - Ek) ;

    double I_VDCC = 0.078*qca_Nav*pow(m_V,2)*h_V*g_Lvdcc*V_spine*(Ca_Cytosol - Ca_Extracellular*exp(-0.078*V_spine)/(1 - exp(-0.078*V_spine))) ;
    double J_Extrusion = (Size_M/Size_C)*P*(Ca_Cytosol - Ca_threshold)*(Ca_Cytosol>Ca_threshold) ;
    double J_Extrusion_D = (1/27.98)*P*(Ca_D_Cytosol - Ca_threshold)*(Ca_D_Cytosol>Ca_threshold) ;
    double J_VDCC = 5182*I_VDCC;

    //double eta_Ca = 1/(Pw4 + Pw1/(Pw2 + pow((Ca_Cytosol-0.05),Pw3))) ;
    plasticity_rate = (-weight + Theta_Ca)/tau_w ;
    //Storing all functional values


    engine::neuron_value(index,"J_IP3R_flux",J_IP3R_flux);
    engine::neuron_value(index,"Flux_NMDA",Flux_NMDA);
    engine::neuron_value(index,"I_NMDA",I_NMDA);
    engine::neuron_value(index,"Eca",Eca);
    engine::neuron_value(index,"I_AMPA",I_AMPA);
    engine::neuron_value(index,"I_SK",I_SK);
    engine::neuron_value(index,"Calcium_flux_total",(/* - J_PACa*/ - J_Ca_d /*- J_CGbinding*/ - J_CaD28k_med - J_CD28k_high - (KFlux_ERM_Cytosol * J_ER_leak_flux) - (KFlux_ERM_Cytosol * J_SERCA_flux) + (KFlux_PM_Cytosol * J_flux0) - (KFlux_ERM_Cytosol * J_IP3R_flux)));


    // ODE
    dxdt[Gbeta_M_index] = ( - J_G1beta + J_NE_G - J_G2beta - J_reconstitution);    // rate for Gbeta_M
    dxdt[GaGTP_PLC_M_index] = (J_NE_GaP + J_PLCassoc -J_GTPase_GaP);    // rate for GaGTP_PLC_M
    //dxdt[SK_M_index] =  0;    // rate for

    dxdt[RL_M_index] = ( - J_G2beta + J_L1 - J_G2);    // rate for RL_M
    dxdt[R_M_index] = ( - J_G1beta - J_G1 - J_L1);    // rate for R_M
    dxdt[RLG_GDP_M_index] = ( - J_NE_RLG + J_L2 + J_G2);    // rate for RLG_GDP_M
    dxdt[RLGbeta_M_index] = (J_L2beta + J_G2beta + J_NE_RLG);    // rate for RLGbeta_M
    dxdt[IP3_PH_YFP_C_index] = J_PH_YFP_IP3;    // rate for IP3_PH_YFP_C
    dxdt[oxoM_EX_index] = ( - (UnitFactor_uM_um3_molecules_neg_1 * KFlux_M_EX * J_L2beta) + J_Oxo_appl - (UnitFactor_uM_um3_molecules_neg_1 * KFlux_M_EX * J_L1) - (UnitFactor_uM_um3_molecules_neg_1 * KFlux_M_EX * J_L2));    // rate for oxoM_EX
    dxdt[PH_CFP_PIP2_M_index] = J_PH_CFP_PIP2;    // rate for PH_CFP_PIP2_M
    dxdt[PH_CFP_C_index] = ( - (UnitFactor_uM_um3_molecules_neg_1 * KFlux_M_C * J_PH_CFP_PIP2) - J_PH_CFP_IP3);    // rate for PH_CFP_C
    dxdt[PI4P_M_index] = ( - J_PIP5K_5Pase + J_PI4K_4Pase - J_PLC_on_PI4P);    // rate for PI4P_M
    dxdt[RGbeta_M_index] = (J_G1beta - J_L2beta + J_NE_RG);    // rate for RGbeta_M
    dxdt[Ga_GDP_M_index] = (J_PLCdiss - J_reconstitution + J_GTPase_Ga);    // rate for Ga_GDP_M
    dxdt[PH_YFP_PIP2_M_index] = J_PH_YFP_PIP2;    // rate for PH_YFP_PIP2_M
    dxdt[DAG_M_index] = (J_PIP2hydr + J_PLC_on_PI4P);    // rate for DAG_M
    dxdt[IP3_C_index] = ( - J_PH_YFP_IP3 - J_IP3Pase + (UnitFactor_uM_um3_molecules_neg_1 * KFlux_M_C * J_PIP2hydr) - J_PH_CFP_IP3);    // rate for IP3_C
    dxdt[GaGTP_M_index] = (J_NE_G - J_PLCassoc + J_NE_RLG + J_NE_RG - J_GTPase_Ga);    // rate for GaGTP_M
    dxdt[PIP2_M_index] = ( - J_PH_CFP_PIP2 - J_PH_YFP_PIP2 + J_PIP5K_5Pase /*- J_PIP2bindSK*/ - J_PIP2hydr);    // rate for PIP2_M
    
    dxdt[D28kB_high_Cytosol_index] =        (J_CD28k_high /*- J_D28kB_high_d*/);    // rate for D28kB_high_Cytosol
    dxdt[IP3_D_Cytosol_index] =         ( - J_IP3deg - J_IP3degrade_IP3_degr1);    // rate for IP3_D_Cytosol
    dxdt[h_ERM_index] =          - J_reaction0;    // rate for h_ERM
    //dxdt[CGB_Cytosol_index] =       (J_CGbinding /*- J_CGB_d*/);    // rate for CGB_Cytosol
    dxdt[h_D_ERM_index] =        - J_reaction1;    // rate for h_D_ERM
    dxdt[Ca_D_Cytosol_index] =      (-J_Extrusion_D - J_CD28k_highDbinding - /*J_PA_Dbinding -*/ J_D28kBDbinding -/* J_CG_Dbinding -*/ J_Ca_deg - (KFlux_ERM_Cytosol * J_SERCA_fluxD) - (KFlux_ERM_Cytosol * J_IP3R_fluxD) + (KFlux_PM_Cytosol * J_flux1) - (KFlux_ERM_Cytosol * J_ER_leak_fluxD));    // rate for Ca_D_Cytosol
    dxdt[D28kB_Cytosol_index] =         ( - J_D28kB_d + J_CaD28k_med);    // rate for D28kB_Cytosol
    //dxdt[PABMg_D_Cytosol_index] =       (J_PAMg_PA_MgD - J_PABMg_deg);    // rate for PABMg_D_Cytosol
    //dxdt[IP3_Cytosol_index] =       ( - J_IP3_d - J_IP3degrade_IP3_degr + J_pulses);    // rate for IP3_Cytosol
    //dxdt[PA_D_Cytosol_index] =      ( - J_PA_Dbinding - J_PAMg_PA_MgD - J_PA_deg);    // rate for PA_D_Cytosol
    dxdt[D28kB_D_Cytosol_index] =       ( - J_D28kB_deg + J_D28kBDbinding);    // rate for D28kB_D_Cytosol
    //dxdt[CG_D_Cytosol_index] =      ( - J_CG_deg - J_CG_Dbinding);    // rate for CG_D_Cytosol
    dxdt[D28kB_high_D_Cytosol_index] =      (J_CD28k_highDbinding - J_D28kB_high_deg);    // rate for D28kB_high_D_Cytosol
    //dxdt[PABCa_Cytosol_index] =         (J_PACa /*- J_PABCa_d*/);    // rate for PABCa_Cytosol
    dxdt[D28k_high_Cytosol_index] =         ( /*- J_D28k_high_d*/ - J_CD28k_high);    // rate for D28k_high_Cytosol
    dxdt[D28k_high_D_Cytosol_index] =       ( - J_CD28k_highDbinding - J_D28k_high_deg);    // rate for D28k_high_D_Cytosol
    //dxdt[CG_Cytosol_index] =        ( /*- J_CG_d */- J_CGbinding);    // rate for CG_Cytosol
    //dxdt[PABMg_Cytosol_index] =         ( /*- J_PABMg_d */+ J_PAMg_PA_Mg);    // rate for PABMg_Cytosol
    //dxdt[CGB_D_Cytosol_index] =         ( - J_CGB_deg + J_CG_Dbinding);    // rate for CGB_D_Cytosol
    dxdt[D28k_Cytosol_index] =      ( - J_CaD28k_med /*- J_D28k_d*/);    // rate for D28k_Cytosol
    //dxdt[PABCa_D_Cytosol_index] =       (J_PA_Dbinding - J_PABCa_deg);    // rate for PABCa_D_Cytosol
    dxdt[Ca_Cytosol_index] =        ( -J_Extrusion -J_VDCC/Size_Cytosol /* - J_PACa*/ - J_Ca_d /*- J_CGbinding - J_CaD28k_med - J_CD28k_high*/ - Flux_NMDA/Size_Cytosol - (KFlux_ERM_Cytosol * J_ER_leak_flux) - (KFlux_ERM_Cytosol * J_SERCA_flux) + (KFlux_PM_Cytosol * J_flux0) - (KFlux_ERM_Cytosol * J_IP3R_flux));    // rate for Ca_Cytosol
    //dxdt[PA_Cytosol_index] =        ( - J_PACa /*- J_PA_d*/ - J_PAMg_PA_Mg);    // rate for PA_Cytosol
    dxdt[V_spine_index] = (-I_dendrite -I_NMDA -I_AMPA -I_SK -I_VDCC)/C_spine;
    dxdt[D28k_D_Cytosol_index] =        ( - J_D28k_deg - J_D28kBDbinding);    // rate for D28k_D_Cytosol

    dxdt[s_index] = (s_sat - s)/tau_s;  // rate of SK activation
    dxdt[weight_index] = plasticity_rate ; //rate os synaptic weight change

    dxdt[m_V_index] = (m_inf - m_V)/tau_m;
    dxdt[h_V_index] = (h_inf - h_V)/tau_h;

    //dxdt[Lpo_index] = +(1-H_SK)+-H_SK/tau_Hpo;

  }
};

} // insilico

#endif
