--  
-- Protocol for Application  Database Interface (PADI)
--
-- Copyright 1995, 1996  Bank of Canada.
--
-- The user of this software has the right to use, reproduce and distribute it.
-- Bank of Canada makes no warranties with respect to the software or its 
-- fitness for any particular purpose. The software is distributed by the Bank
-- of Canada solely on an "as is" basis. By using the software, user agrees to 
-- accept the entire risk of using this software.
--

procedure d_sys
--
-- create test sysdb
--
case 1 to *
series g_compare_formula_db_list_str:string by case = "etsusa", "etscdnfor"
series g_etsbcpi_formula_db_list_str:string by case = "etsusa", "etscdnfor"
series g_etsbis_formula_db_list_str:string by case  = "etsusa", "etscdnfor"
series g_etsboptrade_formula_db_list_str:string by case = "etsusa", "etscdnfor"
series g_etscdnfor_formula_db_list_str:string by case   = "etsusa", "etscdnfor"
series g_etsconsumption_formula_db_list_str:string by case = "etsusa", "etscdnfor"
series g_etscost_formula_db_list_str:string by case    = "etsusa", "etscdnfor"
series g_etscpi_formula_db_list_str:string by case     = "etsusa", "etscdnfor"
series g_etsdrixref_formula_db_list_str:string by case = "etsusa", "etscdnfor"
series g_etsfinance_formula_db_list_str:string by case = "etsusa", "etscdnfor"
series g_etsfinstat_formula_db_list_str:string by case = "etsusa", "etscdnfor"
series g_etsindustrial_formula_db_list_str:string by case = "etsusa", "etscdnfor"
series g_etsintdev_formula_db_list_str:string by case     = "etsusa", "etscdnfor"
series g_etsinvestment_formula_db_list_str:string by case = "etsusa", "etscdnfor"
series g_etsippi_formula_db_list_str:string by case       = "etsusa", "etscdnfor"
series g_etslabour_formula_db_list_str:string by case     = "etsusa", "etscdnfor"
series g_etsmfacansim_formula_db_list_str:string by case  = "etsusa", "etscdnfor"
series g_etsmonmrk_formula_db_list_str:string by case     = "etsusa", "etscdnfor"
series g_etsmrtsxref_formula_db_list_str:string by case   = "etsusa", "etscdnfor"
series g_etsmsoi_formula_db_list_str:string by case       = "etsusa", "etscdnfor"
series g_etspwinv_formula_db_list_str:string by case      = "etsusa", "etscdnfor"
series g_etsregiontemp_formula_db_list_str:string by case = "etsusa", "etscdnfor"
series g_etsusa_formula_db_list_str:string by case        = "etsusa", "etscdnfor"
series g_etswealth_formula_db_list_str:string by case     = "etsusa", "etscdnfor"
series g_etswfs_formula_db_list_str:string by case        = "etsusa", "etscdnfor"

END PROCEDURE
