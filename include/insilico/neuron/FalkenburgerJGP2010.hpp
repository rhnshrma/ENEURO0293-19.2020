/*
 neuron/N_FA.hpp - KCNQ regulation with M1 activation (FalkenburgerJGP2010)

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



namespace insilico {

        


class N_FA: public Neuron{
 public:
  void ode_set(state_type &variables, state_type &dxdt, const double t,const unsigned index) {


    unsigned Gbeta_M_index = engine::neuron_index(index,"Gbeta_M");
    unsigned GaGTP_PLC_M_index = engine::neuron_index(index,"GaGTP_PLC_M");
    unsigned KCNQ_M_index = engine::neuron_index(index,"KCNQ_M");
    unsigned RL_M_index = engine::neuron_index(index,"RL_M");
    unsigned R_M_index = engine::neuron_index(index,"R_M");
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
    double Gbeta_M = variables[Gbeta_M_index];
    double GaGTP_PLC_M = variables[GaGTP_PLC_M_index];
    double KCNQ_M = variables[KCNQ_M_index];
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

    // Load Parameters
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
    double mlabfix_F_nmol_ = engine::neuron_value(index,"mlabfix_F_nmol_");
    double mlabfix_T_ = engine::neuron_value(index,"mlabfix_T_");
    double netValence_NE_RG = engine::neuron_value(index,"netValence_NE_RG");
    double K_millivolts_per_volt = engine::neuron_value(index,"K_millivolts_per_volt");
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
    double mlabfix_PI_ = engine::neuron_value(index,"mlabfix_PI_");
    double PI4P_M_init_molecules_um_2 = engine::neuron_value(index,"PI4P_M_init_molecules_um_2");
    double mlabfix_F_ = engine::neuron_value(index,"mlabfix_F_");
    double netValence_NE_RLG = engine::neuron_value(index,"netValence_NE_RLG");
    double RLG_GDP_M_init_molecules_um_2 = engine::neuron_value(index,"RLG_GDP_M_init_molecules_um_2");
    double mlabfix_R_ = engine::neuron_value(index,"mlabfix_R_");
    double K_IP3ase = engine::neuron_value(index,"K_IP3ase");
    double Kf_GTPase_Ga = engine::neuron_value(index,"Kf_GTPase_Ga");
    double Kf_PLCdiss = engine::neuron_value(index,"Kf_PLCdiss");
    double G_GDP_M_init_molecules_um_2 = engine::neuron_value(index,"G_GDP_M_init_molecules_um_2");
    double Kr_NE_G = engine::neuron_value(index,"Kr_NE_G");
    double mlabfix_K_GHK_ = engine::neuron_value(index,"mlabfix_K_GHK_");
    double surface = engine::neuron_value(index,"surface");
    double speed_PIP2_KCNQ = engine::neuron_value(index,"speed_PIP2_KCNQ");
    double mlabfix_N_pmol_ = engine::neuron_value(index,"mlabfix_N_pmol_");
    double R_M_init_molecules_um_2 = engine::neuron_value(index,"R_M_init_molecules_um_2");
    double KD_PH_IP3 = engine::neuron_value(index,"KD_PH_IP3");
    double Size_M = engine::neuron_value(index,"Size_M");
    double Size_C = engine::neuron_value(index,"Size_C");
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
    double KMOLE = engine::neuron_value(index,"KMOLE");
    double K4_Oxo = engine::neuron_value(index,"K4_Oxo");


    unsigned on1 = 0;
    unsigned on2 = 0;
    unsigned off1 = 0;
    unsigned off2 = 0;

    if (t>950){
      on1 = 1;
      off1 = 0;
    }
    else{
      on1 = 0;
      off1 = 1;
    }
    if (t<990){
      on2 = 1;
      off2 = 0;
    }
    else{
      on2 = 0;
      off2 = 1;
    }

    // Functions for calculations
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
    double KCNQ_PIP2_M = ((K_KCNQ_PIP2_M_total - (Size_M * UnitFactor_uM_um3_molecules_neg_1 * KCNQ_M)) / (Size_M * UnitFactor_uM_um3_molecules_neg_1));
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
    double Kf_PIP2bindKCNQ = speed_PIP2_KCNQ;
    double KD_PIP2_KCNQ = pow(KA_PIP2_KCNQ, Hill_binding);
    double Kr_PIP2bindKCNQ = (KD_PIP2_KCNQ * Kf_PIP2bindKCNQ);
    double J_PIP2bindKCNQ = ((Kf_PIP2bindKCNQ * pow(PIP2_M , Hill_binding) * KCNQ_M) - (Kr_PIP2bindKCNQ * KCNQ_PIP2_M));
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

    

    //Storing all functional values
    engine::neuron_value(index,"UnitFactor_uM_um3_molecules_neg_1",UnitFactor_uM_um3_molecules_neg_1);
    engine::neuron_value(index,"K_OxoM_conc_EX_total",K_OxoM_conc_EX_total);
    engine::neuron_value(index,"Kf_PH_YFP_IP3",Kf_PH_YFP_IP3);
    engine::neuron_value(index,"Kr_PH_YFP_IP3",Kr_PH_YFP_IP3);
    engine::neuron_value(index,"K_PH_YFP_C_total",K_PH_YFP_C_total);
    engine::neuron_value(index,"PH_YFP_C",PH_YFP_C);
    engine::neuron_value(index,"J_PH_YFP_IP3",J_PH_YFP_IP3);
    engine::neuron_value(index,"KfL1",KfL1);
    engine::neuron_value(index,"KfL2",KfL2);
    engine::neuron_value(index,"K_EX",K_EX);
    engine::neuron_value(index,"Kr_PH_YFP_PIP2",Kr_PH_YFP_PIP2);
    engine::neuron_value(index,"Kf_PH_YFP_PIP2",Kf_PH_YFP_PIP2);
    engine::neuron_value(index,"J_PH_YFP_PIP2",J_PH_YFP_PIP2);
    engine::neuron_value(index,"K_KCNQ_PIP2_M_total",K_KCNQ_PIP2_M_total);
    engine::neuron_value(index,"KCNQ_PIP2_M",KCNQ_PIP2_M);
    engine::neuron_value(index,"K_GaGDP_PLC_M_total",K_GaGDP_PLC_M_total);
    engine::neuron_value(index,"GaGDP_PLC_M",GaGDP_PLC_M);
    engine::neuron_value(index,"J_NE_GaP",J_NE_GaP);
    engine::neuron_value(index,"K_IP3_PH_CFP_C_total",K_IP3_PH_CFP_C_total);
    engine::neuron_value(index,"IP3_PH_CFP_C",IP3_PH_CFP_C);
    engine::neuron_value(index,"allPH_cytoplasm",allPH_cytoplasm);
    engine::neuron_value(index,"K_G_GDP_M_total",K_G_GDP_M_total);
    engine::neuron_value(index,"G_GDP_M",G_GDP_M);
    engine::neuron_value(index,"J_NE_G",J_NE_G);
    engine::neuron_value(index,"J_PLC_on_PI4P",J_PLC_on_PI4P);
    engine::neuron_value(index,"Kr_PH_CFP_IP3",Kr_PH_CFP_IP3);
    engine::neuron_value(index,"Kf_PH_CFP_PIP2",Kf_PH_CFP_PIP2);
    engine::neuron_value(index,"K_PLC_M_total",K_PLC_M_total);
    engine::neuron_value(index,"PLC_M",PLC_M);
    engine::neuron_value(index,"J_PLCdiss",J_PLCdiss);
    engine::neuron_value(index,"Kf_PIP2bindKCNQ",Kf_PIP2bindKCNQ);
    engine::neuron_value(index,"KD_PIP2_KCNQ",KD_PIP2_KCNQ);
    engine::neuron_value(index,"Kr_PIP2bindKCNQ",Kr_PIP2bindKCNQ);
    engine::neuron_value(index,"J_PIP2bindKCNQ",J_PIP2bindKCNQ);
    engine::neuron_value(index,"allRL",allRL);
    engine::neuron_value(index,"KG2",KG2);
    engine::neuron_value(index,"KG1",KG1);
    engine::neuron_value(index,"KfG1",KfG1);
    engine::neuron_value(index,"KrG1",KrG1);
    engine::neuron_value(index,"Kr_G1beta",Kr_G1beta);
    engine::neuron_value(index,"J_NE_RLG",J_NE_RLG);
    engine::neuron_value(index,"KL2",KL2);
    engine::neuron_value(index,"KrL2",KrL2);
    engine::neuron_value(index,"Kr_L2beta",Kr_L2beta);
    engine::neuron_value(index,"PI_M",PI_M);
    engine::neuron_value(index,"Kf_L2",Kf_L2);
    engine::neuron_value(index,"Kf_L1",Kf_L1);
    engine::neuron_value(index,"Kr_G2",Kr_G2);
    engine::neuron_value(index,"Kr_G1",Kr_G1);
    engine::neuron_value(index,"Kf_Oxo_appl",Kf_Oxo_appl);
    engine::neuron_value(index,"K_RG_GDP_M_total",K_RG_GDP_M_total);
    engine::neuron_value(index,"RG_GDP_M",RG_GDP_M);
    engine::neuron_value(index,"Kr_L2",Kr_L2);
    engine::neuron_value(index,"J_L2",J_L2);
    engine::neuron_value(index,"Kr_L1",Kr_L1);
    engine::neuron_value(index,"J_L1",J_L1);
    engine::neuron_value(index,"Kf_PH_CFP_IP3",Kf_PH_CFP_IP3);
    engine::neuron_value(index,"Kf_G2beta",Kf_G2beta);
    engine::neuron_value(index,"Kr_PH_CFP_PIP2",Kr_PH_CFP_PIP2);
    engine::neuron_value(index,"J_PH_CFP_PIP2",J_PH_CFP_PIP2);
    engine::neuron_value(index,"Kr_G2beta",Kr_G2beta);
    engine::neuron_value(index,"J_G2beta",J_G2beta);
    engine::neuron_value(index,"Kf_PI4K_4Pase",Kf_PI4K_4Pase);
    engine::neuron_value(index,"J_reconstitution",J_reconstitution);
    engine::neuron_value(index,"J_PH_CFP_IP3",J_PH_CFP_IP3);
    engine::neuron_value(index,"J_GTPase_GaP",J_GTPase_GaP);
    engine::neuron_value(index,"allGaPLC",allGaPLC);
    engine::neuron_value(index,"Kr_Oxo_appl",Kr_Oxo_appl);
    engine::neuron_value(index,"J_PI4K_4Pase",J_PI4K_4Pase);
    engine::neuron_value(index,"J_NE_RG",J_NE_RG);
    engine::neuron_value(index,"J_IP3Pase",J_IP3Pase);
    engine::neuron_value(index,"Kf_PIP5K_5Pase",Kf_PIP5K_5Pase);
    engine::neuron_value(index,"Kf_G1beta",Kf_G1beta);
    engine::neuron_value(index,"Kf_L2beta",Kf_L2beta);
    engine::neuron_value(index,"J_PLCassoc",J_PLCassoc);
    engine::neuron_value(index,"Kf_G2",Kf_G2);
    engine::neuron_value(index,"Kf_G1",Kf_G1);
    engine::neuron_value(index,"KFlux_M_EX",KFlux_M_EX);
    engine::neuron_value(index,"J_G1beta",J_G1beta);
    engine::neuron_value(index,"J_L2beta",J_L2beta);
    engine::neuron_value(index,"J_PIP2hydr",J_PIP2hydr);
    engine::neuron_value(index,"J_G2",J_G2);
    engine::neuron_value(index,"J_G1",J_G1);
    engine::neuron_value(index,"J_GTPase_Ga",J_GTPase_Ga);
    engine::neuron_value(index,"OxoM_conc_EX",OxoM_conc_EX);
    engine::neuron_value(index,"J_Oxo_appl",J_Oxo_appl);
    engine::neuron_value(index,"allPIP2M",allPIP2M);
    engine::neuron_value(index,"IM",IM);
    engine::neuron_value(index,"KFlux_M_C",KFlux_M_C);
    engine::neuron_value(index,"J_PIP5K_5Pase",J_PIP5K_5Pase);
    engine::neuron_value(index,"K_C",K_C);
    engine::neuron_value(index,"allRGbeta",allRGbeta);
    engine::neuron_value(index,"allGabg",allGabg);
    engine::neuron_value(index,"Kf_Oxo_appl",Kf_Oxo_appl);
    engine::neuron_value(index,"J_PIP2bindKCNQ",J_PIP2bindKCNQ);




    // ODE
    dxdt[Gbeta_M_index] = ( - J_G1beta + J_NE_G - J_G2beta - J_reconstitution);    // rate for Gbeta_M
    dxdt[GaGTP_PLC_M_index] = (J_NE_GaP + J_PLCassoc - J_GTPase_GaP);    // rate for GaGTP_PLC_M
    dxdt[KCNQ_M_index] =  - J_PIP2bindKCNQ;    // rate for KCNQ_M
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
    dxdt[PIP2_M_index] = ( - J_PH_CFP_PIP2 - J_PH_YFP_PIP2 + J_PIP5K_5Pase - J_PIP2bindKCNQ - J_PIP2hydr);    // rate for PIP2_M
  }
};

} // insilico

#endif
