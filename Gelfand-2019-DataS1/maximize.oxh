/*--------------------------------------------------------------------------
 * maximize.h - definitions/declarations for numerical differentiation
 *              and function maximization
 *
 *       (C) Jurgen Doornik 1994-2001
 *
 *--------------------------------------------------------------------------*/

#ifndef MAXIMIZE_INCLUDED
#define MAXIMIZE_INCLUDED

enum
{
    MAX_CONV, MAX_WEAK_CONV, MAX_MAXIT, MAX_LINE_FAIL, MAX_FUNC_FAIL, MAX_NOCONV
};
MaxMonitor(const sMethod, const cIter, const vP, const vScore, const vDelta,
    const dFuncInit, const dFunc, const steplen, const condhes, const ret_val,
	const bCompact);
MaxGetParams(const args);
MaxGetReturn(const iResult, const cIter, const objMaxctrl);

MaxBFGS(const func, const avP, const adFunc, const amHessian, const fNumDer, ...);
MaxNewton(const func, const avP, const adFunc, const amHessian, const fNumDer, ...);
MaxSimplex(const func, const avP, const adFunc, vDelta, ...);
MaxConvergenceMsg(const iCode);
MaxControlEps(const dEps1, const dEps2);
MaxControl(const mxIter, const iPrint, ...);
GetMaxControl();
GetMaxControlEps();

Num1Derivative(const func, vP, const avScore);
Num1Derivative_parallel(const func, vP, const avScore);
Num2Derivative(const func, vP, const amHessian);
Num2Derivative_parallel(const func, vP, const amHessian);
NumJacobian(const func, vU, const amJacobian);
NumJacobianEx(const func, vU, const vF, const fCentral, const amJacobian);
NumJacobianMul(const func, const vU, const vF, const mD, const fCentral,
	const amJacobian);

class CMaxControl
{
	CMaxControl(const iOptions = 0);
	SetOptions(const iOptions);
	SetEps(const dEps1, const dEps2=-1);
	GetEps();
	SetControl(const mxIter, const iPrint=-1, const bCompact=-1);
	GetControl();
	GetResult();
	GetIterationCount();
	SetResult(const iResult);
	SetIterationCount(const cIter);

	decl m_mxIter;							    /* maximum no of iterations */
	decl m_dEps1, m_dEps2;  			            /* convergence criteria */
	decl m_iPrint;    		   /* print results every s_iPrint'th iteration */
	decl m_bCompact;	 	    /* TRUE: print iterations in compact format */
	decl m_fnNum1der;					 		/* numerical score function */

	decl m_iResult;								 	  /* convergence result */
	decl m_cIter;									 	 /* iteration count */
public:
	enum
	{	PARALLEL_SCORE = 1
	};
}


#endif /* MAXIMIZE_INCLUDED */
