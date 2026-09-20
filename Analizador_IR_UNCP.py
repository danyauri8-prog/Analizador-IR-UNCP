# ANALIZADOR VECTORIAL IR - UNCP
# Cristopher Daniel Yauri Navarro
import math, numpy as np, matplotlib.pyplot as plt
C=299792458.0; H=6.62607015e-34; NA=6.02214076e23
def calcular_proyeccion(dm, eh):
    dm=np.array(dm,float); eh=np.array(eh,float)
    eh=eh/np.linalg.norm(eh)
    prod=np.dot(dm,eh); mag=np.linalg.norm(dm)
    ang=math.degrees(math.acos(np.clip(prod/mag,-1,1))) if mag!=0 else 0
    return prod,mag,ang,eh
def convertir(wn):
    wm=wn*100; f=C*wm; lam=1/wm; Ej=H*f; EkJ=(Ej*NA)/1000
    return f/1e12, lam*1e6, EkJ
print("EJEMPLO 1715 cm-1 = 20.516 kJ/mol")
dm=np.array([0.8,0.5,0.3])*1e-30; ed=np.array([1,0.2,0.1])
prod,mag,ang,en=calcular_proyeccion(dm,ed)
print(f"Delta magnitude {mag:.3e} C·m Angulo {ang:.2f}°")
print(f"Delta·e = {prod:.3e} C·m")
wn=1715; fT, lamU, EkJ=convertir(wn)
print(f"Para {wn} cm-1: {fT:.3f} THz, {lamU:.3f} um, {EkJ:.3f} kJ/mol")
x=np.linspace(4000,400,1000)
g=lambda x,mu,s,a: a*np.exp(-0.5*((x-mu)/s)**2)
y=g(x,1715,20,1)+g(x,2950,30,0.6)+g(x,3300,40,0.4)
plt.plot(x,y); plt.gca().invert_xaxis()
plt.title("Espectro IR - UNCP"); plt.xlabel("cm-1"); plt.ylabel("Abs"); plt.grid(True,alpha=0.3); plt.show()
