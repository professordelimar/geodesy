"""Funções básicas para sistemas de tempo usados em Geodésia."""
import numpy as np
from datetime import datetime, timedelta
def julian_date(year,month,day,hour=0,minute=0,second=0.0):
    y=int(year); m=int(month); d=float(day)+(hour+minute/60+second/3600)/24
    if m<=2: y-=1; m+=12
    A=int(y/100); B=2-A+int(A/4); return int(365.25*(y+4716))+int(30.6001*(m+1))+d+B-1524.5
def modified_julian_date(year,month,day,hour=0,minute=0,second=0.0): return julian_date(year,month,day,hour,minute,second)-2400000.5
def datetime_to_julian_date(dt): return julian_date(dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second+dt.microsecond/1e6)
def calendar_to_decimal_year(year,month,day,hour=0,minute=0,second=0.0):
    dt=datetime(year,month,day,hour,minute,int(second)); start=datetime(year,1,1); end=datetime(year+1,1,1); return year+(dt-start).total_seconds()/(end-start).total_seconds()
def decimal_year_to_datetime(decimal_year):
    year=int(np.floor(decimal_year)); frac=decimal_year-year; start=datetime(year,1,1); end=datetime(year+1,1,1); return start+timedelta(seconds=frac*(end-start).total_seconds())
def gps_week_seconds(dt):
    ep=datetime(1980,1,6); sec=(dt-ep).total_seconds(); week=int(sec//(7*86400)); return week,sec-week*7*86400
def greenwich_sidereal_time_approx(dt):
    JD=datetime_to_julian_date(dt); T=(JD-2451545.0)/36525.0; theta=280.46061837+360.98564736629*(JD-2451545.0)+0.000387933*T*T-T**3/38710000.0; return theta%360.0
