

#ifndef ReportMCMCH_INCLUDED
#define ReportMCMCH_INCLUDED

class ReportMCMC
{

	decl m_sversion; // Version
	decl m_cdim, m_cRep;	 // dimension of samples, # of samples
	decl m_msample, m_sname; // samples, variable name of sample
	decl m_dBm, m_sfmt, m_sfigtype;		 // bandwidth, output format, figure type
	decl m_asname, m_soutfilename; // name of variables, name of outputfile
	decl m_fcorr; // flag: 1: Report Correlation Matrix for MCMC output, 0: No report
	decl m_cmaxnumfig, m_cpage; // control output of figures
	ReportMCMC(mdata);// constructor

	Report();
	
	DrawDensityPlot();
	DrawACFPlot();
	DrawSamplePath();
	CalculateStatistics();
	OutputFile();

	SetBandwidth(const dBm);	// change the bandwidth
	SetFormat(const sformat);	// change the output format
	SetFigureType(const sfigtype); // change the figure type
	SetOutfileName(const soutfilename); // change outputfile name
	SetVarNames(const asname); // Set variable names
	SetCorrOutput(const fcorr);	// suppress correlation output when fcorr = 0;
	SetMaxNumFig(const cmaxnumfig);	// Set Maximum Number of Figures per page 
	
	fTsvar(const mX, const dBm);
	fTsvar_Batch(const mX, const dBm);
	fGeweke(const mX, const dBm);
	fCI(const mX);
	/*
	 fTsvar calculates a variance of time series.
	 fGeweke calculates p-value for Convergence (Geweke's Method)
	 H0: Convergence, H1: No Convergence
	*/

	
};// the end of class definition

#endif
