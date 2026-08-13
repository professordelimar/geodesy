"""Colocação por mínimos quadrados em versão inicial."""
import numpy as np
def empirical_covariance(values,max_lag=None):
    v=np.asarray(values,dtype=float)-np.nanmean(values); n=len(v); max_lag=n//4 if max_lag is None else max_lag; return np.array([np.nanmean(v[:n-lag]*v[lag:]) for lag in range(max_lag+1)])
def covariance_model_gaussian(distance,sigma2=1.0,correlation_length=1.0): return sigma2*np.exp(-(np.asarray(distance,dtype=float)/correlation_length)**2)
def covariance_model_exponential(distance,sigma2=1.0,correlation_length=1.0): return sigma2*np.exp(-np.asarray(distance,dtype=float)/correlation_length)
def distance_matrix(x1,y1,x2=None,y2=None):
    x1=np.asarray(x1,dtype=float).ravel(); y1=np.asarray(y1,dtype=float).ravel(); x2=x1 if x2 is None else np.asarray(x2,dtype=float).ravel(); y2=y1 if y2 is None else np.asarray(y2,dtype=float).ravel(); return np.sqrt((x1[:,None]-x2[None,:])**2+(y1[:,None]-y2[None,:])**2)
def covariance_matrix(x,y,model='gaussian',sigma2=1.0,correlation_length=1.0,noise=0.0):
    D=distance_matrix(x,y); C=covariance_model_gaussian(D,sigma2,correlation_length) if model=='gaussian' else covariance_model_exponential(D,sigma2,correlation_length); return C+noise*np.eye(C.shape[0])
def least_squares_collocation(x_obs,y_obs,values_obs,x_pred,y_pred,sigma2=1.0,correlation_length=1.0,noise=1e-6,model='gaussian'):
    x_obs=np.asarray(x_obs,dtype=float).ravel(); y_obs=np.asarray(y_obs,dtype=float).ravel(); values_obs=np.asarray(values_obs,dtype=float).ravel(); x_pred=np.asarray(x_pred,dtype=float).ravel(); y_pred=np.asarray(y_pred,dtype=float).ravel()
    D_oo=distance_matrix(x_obs,y_obs); D_po=distance_matrix(x_pred,y_pred,x_obs,y_obs)
    if model=='gaussian': C_oo=covariance_model_gaussian(D_oo,sigma2,correlation_length); C_po=covariance_model_gaussian(D_po,sigma2,correlation_length)
    else: C_oo=covariance_model_exponential(D_oo,sigma2,correlation_length); C_po=covariance_model_exponential(D_po,sigma2,correlation_length)
    return C_po@np.linalg.solve(C_oo+noise*np.eye(len(x_obs)),values_obs)
def collocation_error_variance(x_obs,y_obs,x_pred,y_pred,sigma2=1.0,correlation_length=1.0,noise=1e-6,model='gaussian'):
    x_obs=np.asarray(x_obs,dtype=float).ravel(); y_obs=np.asarray(y_obs,dtype=float).ravel(); x_pred=np.asarray(x_pred,dtype=float).ravel(); y_pred=np.asarray(y_pred,dtype=float).ravel(); D_oo=distance_matrix(x_obs,y_obs); D_po=distance_matrix(x_pred,y_pred,x_obs,y_obs)
    if model=='gaussian': C_oo=covariance_model_gaussian(D_oo,sigma2,correlation_length); C_po=covariance_model_gaussian(D_po,sigma2,correlation_length)
    else: C_oo=covariance_model_exponential(D_oo,sigma2,correlation_length); C_po=covariance_model_exponential(D_po,sigma2,correlation_length)
    C_oo+=noise*np.eye(len(x_obs)); return sigma2-np.sum(C_po*np.linalg.solve(C_oo,C_po.T).T,axis=1)
