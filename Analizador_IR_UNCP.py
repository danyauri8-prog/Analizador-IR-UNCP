# ANALIZADOR IR - UNCP - Yauri Navarro
import math, numpy as np, matplotlib.pyplot as plt
C=299792458.0; H=6.62607015e-34; NA=6.02214076e23
dm=np.array([0.8,0.5,0.3])*1e-30
ed=np.array([1,0.2,0.1])
ed=ed/np.linalg.norm(ed)
prod=np.dot(dm,ed); mag=np.linalg.norm(dm)
ang=math.degrees(math.acos(prod/mag))
print(f"EJEMPLO 1715 cm-1 = 20.516 kJ/mol")
print(f"Delta magnitude {mag:.3e} C·m")
print(f"Angulo {ang:.2f}° Delta·e {prod:.3e}")
wn=1715; wm=wn*100; f=C*wm; lam=1/wm; E=(H*f*NA)/1000
print(f"{wn} cm-1 -> {f/1e12:.2f} THz, {lam*1e6:.3f} um, {E:.3f} kJ/mol")
x=np.linspace(4000,400,1000)
g=lambda x,mu,s,a: a*np.exp(-0.5*((x-mu)/s)**2)
y=g(x,1715,20,1)+g(x,2950,30,0.6)+g(x,3300,40,0.4)
plt.plot(x,y); plt.gca().invert_xaxis()
plt.title("Espectro IR UNCP"); plt.xlabel("cm-1"); plt.grid(True,alpha=0.3); plt.show()
