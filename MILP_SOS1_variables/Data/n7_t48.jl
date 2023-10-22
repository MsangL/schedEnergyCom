T=48
Pc=[3.7;5.000;3.700;0.0;0.0;0.0;0.0]                   # Puissances kW
Pd=[3.7;5.000;3.700;0.0;0.0;0.0;0.0]                   # Puissances kW
M=2                                                  # maximum de schedule par tâche
fi=10                                                # nombre max de cycle des battéries
P_sous=[6 36 6 9 9 6 9 ]                             # puissance souscrite convertir
JA=2
JB=1
delta=1800
gamma=[9.8;9.8;9.8;0.0;0.0;0.0;0.0]                # Kwh
n=7                                                  # The community's size
B=1                                                  # Maximum number of 
dist=0.5                                             # chaque 30 min 1/2h
c=0.975
d=0.975
xi=zeros(n)                      # charge initiale des batteries
eta=0.99                                             # perte automatique
alpha=1                                              # pas de perte suivant la distance
M=2                                                  # maximum de schedule par tâche
fi=10                                                # nombre max de cycle des battéries
L=1
Pedf=0.1685 # prix d'achat au réseau edf euro/kWh
Pac=0.14     # prix d'achat dans la communauté euro/kWh
Pvc=0.12     # Prix de vente dans la communauté euro/kWh
oui=0.065   # prix de vente à planète oui euro/kWh
Pv=0.1        # prix de vente à edf lorsque pas dans une communauté euro/kWh

Nbp=3
u=1
tt=0

JobA= Array{Float64}(undef, (n,JA,Nbp,u))

JobA[:,:,1,1]=[1.0 1.0 ; 0.0 0.0 ; 1.0 1.0 ; 1.0 1.0 ;1.0 1.0 ;0.0 1.0 ;1.0 1.0 ]
JobA[:,:,2,1]=[1.0 0.0 ; 0.0 0.0 ; 0.0 0.0 ; 0.0 0.0 ;1.0 0.0 ;1.0 0.0 ;1.0 0.0 ]
JobA[:,:,3,1]=[1.0 0.0 ; 0.0 0.0 ; 1.0 0.0 ; 0.0 0.0 ;0.0 0.0 ;1.0 0.0 ;1.0 0.0 ]

pow=4       # nombre max de puissance
Pa= Array{Float64}(undef, (n,JA,Nbp,pow))

Pa[:,:,:,4].=3.0
Pa[:,:,:,3].=2.0
Pa[:,:,:,2].=1.0
Pa[:,:,:,1].=0.0

Pb= zeros(n,T,JB,M)

tmin=[19 0 16 17 20 22 19;22 0 19 22 23 22 19;16 0 19 17 20 22 19]
tmax=[22 24 19 20 21 24 21;24 24 24 23 24 24 21;20 21 24 20 22 24 21]

Cm=[2970000 2970000 5940000 3960000 8250000 6105000 4950000;4950000 2970000 2970000 6600000 3300000 2970000 3300000;5940000 2970000 2970000 3960000 3960000 3300000 4620000]
U=[12.0 12.0 24.0 16.0 33.3 24.6 20.0;20.0 12.0 12.0 20.0 13.3 0.0 13.3;24.6 12.0 12.0 16.0 16.0 13.3 18.6]

td=[20 1 20 20 34 18 20;20 1 20 20 34 18 20;20 1 20 20 34 18 20]
tf=[30 48 40 40 47 34 42;30 48 40 40 47 34 42;30 48 40 40 47 34 42]

Yout=[4.4 4.1 3.8 3.55 3.3 3.0 2.7 3.0 3.3 4.1 4.9 4.55 4.2 4.75 5.3 5.35 5.4 6.7 8.0 8.25 8.5 9.25 10.0 10.55 11.1 11.45 11.8 12.1 12.4 12.25 12.1 10.65 9.2 8.6 8.0 7.75 7.5 7.1 6.7 7.05 7.4 8.2 9.0 9.4 9.8 9.75 9.7 9.7]

ti=[15 17 18 19 10 10 8;15 17 18 19 10 10 8;15 17 18 19 10 10 8]
y0=[12 14 16 12 14 12 15;15 14 12 11 14 12 8;12 14 10 12 14 12 10]

sc=28
ec=36
temp_min=55
temp_max=60
Ta=17
S=[1.5 7.0 2.0 3.75 2.4 2.0 2.6;1.5 7.0 2.0 3.75 2.4 2.0 2.6;1.5 7.0 2.0 3.75 2.4 2.0 2.6]
Masse=[75 350 100 200 150 100 150;75 350 100 200 150 100 150;75 350 100 200 150 100 150]
K=1
cpo=1
r=1.2

prod1=[0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.008 0.036 0.144 0.124 0.41 0.342 0.31 0.276 0.042 0.044 0.164 0.008 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.014 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.13 0.296 0.184 0.174 0.126 0.03 0.034 0.01 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.122 0.714 0.896 1.12 1.3 1.32 1.406 1.316 1.21 1.004 0.708 0.112 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
 

p=[0.084 0.0254 0.0318 0.0872 0.1668 0.0994 0.0798 0.0318 0.0308 0.0346 0.0944 0.098 0.1468 0.1308 0.197 0.2942 0.3656 0.319 0.245 0.226 0.1574 0.148 0.1576 0.0872 0.0898 0.172 0.5422 0.5154 0.5958 0.3216 0.1418 0.0868 0.1358 0.1694 0.22 0.228 0.222 0.22 0.2426 0.2384 0.2504 0.2654 0.2194 0.2346 0.0902 0.026 0.0334 0.0668;1.164 1.106 1.56 1.008 1.052 1.524 1.1 1.514 1.106 1.184 1.382 1.07 1.35 1.088 1.676 1.148 1.198 1.504 2.584 3.062 2.232 1.646 2.008 1.4 1.262 0.916 0.98 0.94 1.438 1.208 1.382 1.158 0.878 1.448 0.968 0.952 1.566 0.942 0.878 1.082 0.846 1.306 1.06 0.944 1.476 0.992 0.99 1.16;0.2728 0.2568 0.2244 0.1944 0.2162 0.2388 0.1604 0.149 0.1848 0.1986 0.1588 0.1842 0.194 0.1956 0.313 0.1538 0.2198 0.207 0.1596 0.3614 0.2094 0.1668 0.1612 0.2572 0.1566 0.1598 0.1686 0.1594 0.1572 0.1562 0.16 0.1568 0.405 0.3454 0.2426 0.1272 0.154 0.2306 0.1268 0.1478 0.1908 0.133 0.1448 0.1578 0.1452 0.1422 0.134 0.1588;0.2344 0.1588 0.156 0.2732 0.163 0.1488 0.3042 0.3506 0.3986 0.4288 0.2598 0.2212 0.2234 0.21 0.3334 0.5226 0.4044 0.2076 0.2314 0.1806 0.2912 0.353 0.2654 0.0342 0.0676 0.3264 0.5472 0.3282 0.1158 0.1116 0.3222 0.143 0.117 0.248 0.1568 0.1586 0.3962 0.294 0.2102 0.2068 0.177 0.295 0.3848 0.151 0.146 0.1102 0.2288 0.1032;0.1156 0.1348 0.1846 0.122 0.1388 0.1556 0.108 0.1752 0.1312 0.1242 0.1698 0.135 0.1478 0.1572 0.1664 0.3682 0.327 0.393 0.479 0.4054 0.6136 0.5056 0.4188 0.3746 0.2332 0.177 0.0714 0.0806 0.0738 0.127 0.1436 0.1262 0.1362 0.1662 0.1808 0.153 0.118 0.1734 0.229 0.381 0.4476 0.121 0.152 0.0954 0.0548 0.075 0.0776 0.0838;0.0358 0.0448 0.0382 0.0308 0.0358 0.0304 0.0416 0.0384 0.0448 0.0312 0.0336 0.0372 0.1092 0.0646 0.062 0.1258 0.0868 0.0434 0.035 0.0028 0.0002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.025 0.0276 0.0344 0.0466 0.0654 0.0432 0.0548 0.045 0.0496 0.0534 0.0474 0.0524 0.0454 0.0586 0.0432 0.0462 0.0378 0.0328;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
beta=0.15
omega=[1 1 0 1 1 1 0]
maxIterGC=10


