"""Modelos simples para aplicações geodinâmicas e séries temporais geodésicas."""
import numpy as np
def linear_velocity_model(t,x0,v,t0=0.0): return x0+v*(np.asarray(t,dtype=float)-t0)
def estimate_linear_velocity(t,x):
    t=np.asarray(t,dtype=float); x=np.asarray(x,dtype=float); A=np.column_stack([np.ones_like(t),t]); coef,*_=np.linalg.lstsq(A,x,rcond=None); return coef[0],coef[1]
def seasonal_signal(t,offset=0,trend=0,annual_sin=0,annual_cos=0,semi_sin=0,semi_cos=0):
    t=np.asarray(t,dtype=float); return offset+trend*t+annual_sin*np.sin(2*np.pi*t)+annual_cos*np.cos(2*np.pi*t)+semi_sin*np.sin(4*np.pi*t)+semi_cos*np.cos(4*np.pi*t)
def fit_seasonal_signal(t,x):
    t=np.asarray(t,dtype=float); x=np.asarray(x,dtype=float); A=np.column_stack([np.ones_like(t),t,np.sin(2*np.pi*t),np.cos(2*np.pi*t),np.sin(4*np.pi*t),np.cos(4*np.pi*t)]); coef,*_=np.linalg.lstsq(A,x,rcond=None); return coef
def postseismic_log_decay(t,A,tau,t0=0.0): return A*np.log1p(np.maximum(np.asarray(t,dtype=float)-t0,0.0)/tau)
def sea_level_trend(t,sea_level): return estimate_linear_velocity(t,sea_level)
