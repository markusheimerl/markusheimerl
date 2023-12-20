#/*
#* Markus Heimerl
#* 2021, OTH Regensburg
#* markus (at) markusheimerl (dot) com
#*/

import math
from math import sin, cos, tan, sqrt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from mpl_toolkits.mplot3d import Axes3D

pylab.rcParams['figure.figsize'] = 8, 8

class Drohne3D():
    
    def __init__(self, k_f, k_m, m, L, i_x, i_y, i_z, omega_quadrat_min, omega_quadrat_max):
        self.k_f = k_f
        self.k_m = k_m
        self.m = m
        self.l = (L * sqrt(2)) / 2 # siehe Abbildung "3D Drohne von oben" (Abbildungsverzeichnis)
        self.i_x = i_x
        self.i_y = i_y
        self.i_z = i_z
        self.omega_quadrat_min = omega_quadrat_min
        self.omega_quadrat_max = omega_quadrat_max
        
        # siehe Gleichung 21
        self.X=np.array([
            # x, y, z, phi, theta, psi, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            # x_punkt, y_punkt, z_punkt, p, q, r
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        self.omega = np.array([0.0, 0.0, 0.0, 0.0])
        self.g = 9.81
    
    # Positionen (Im Weltkoordinatensystem)
    @property
    def x(self):
        return self.X[0]

    @property
    def y(self):
        return self.X[1]

    @property
    def z(self):
        return self.X[2]   
    
    # Eulerwinkel [rad] (Im Weltkoordinatensystem)
    @property
    def phi(self):
        return self.X[3]

    @property
    def theta(self):
        return self.X[4]

    @property
    def psi(self):
        return self.X[5]

    # Körperachsengeschwindigkeiten [rad / s]
    @property 
    def p(self):
        return self.X[9]

    @property
    def q(self):
        return self.X[10]

    @property 
    def r(self):
        return self.X[11]
        
    # Kräfte der vier Propeller [N]
    @property
    def f_1(self):
        # siehe Gleichung 2
        return self.k_f*self.omega[0]**2

    @property 
    def f_2(self):
        # siehe Gleichung 2
        return self.k_f*self.omega[1]**2

    @property 
    def f_3(self):
        # siehe Gleichung 2
        return self.k_f*self.omega[2]**2

    @property 
    def f_4(self):
        # siehe Gleichung 2
        return self.k_f*self.omega[3]**2

    @property
    def f_gesamt(self):
        # siehe Gleichung 25
        return self.f_1 + self.f_2 + self.f_3 + self.f_4
        
    # Momente der vier Propeller [N * m]
    @property
    def m_1(self):
        # siehe Gleichung 43
        return self.k_m * self.omega[0]**2
        
    @property
    def m_2(self):
        # siehe Gleichung 43
        return self.k_m * self.omega[1]**2

    @property
    def m_3(self):
        # siehe Gleichung 43
        return self.k_m * self.omega[2]**2

    @property
    def m_4(self):
        # siehe Gleichung 43
        return self.k_m * self.omega[3]**2
        
    # Momente um die Axen [N * m]
    @property
    def m_x(self):
        # siehe Gleichung 42
        return self.l*((self.f_1 + self.f_3) - (self.f_2  + self.f_4))

    @property
    def m_y(self):
        # siehe Gleichung 42
        return self.l*((self.f_1 + self.f_2) - (self.f_3 + self.f_4))

    @property
    def m_z(self):
        # siehe Gleichung 42
        return (self.m_1 + self.m_4) - (self.m_2 + self.m_3)
        
    def setze_rotationsgeschwindigkeiten(self, f_gesamt_ziel, m_x_ziel, m_y_ziel, m_z_ziel):
        # siehe Gleichung 44
        inv_mat = np.array([[0.25, 0.25, 0.25, 0.25], [0.25, -0.25, 0.25, -0.25], [0.25, 0.25, -0.25, -0.25], [0.25, -0.25, -0.25, 0.25]])
        ziel_vek = np.array([(f_gesamt_ziel/self.k_f), (m_x_ziel/(self.l*self.k_f)), (m_y_ziel/(self.l*self.k_f)), (m_z_ziel/(self.k_m))]).T
        omega_ziel_vek = np.matmul(inv_mat, ziel_vek)

        self.omega[0] = np.sqrt(min(self.omega_quadrat_max, max(self.omega_quadrat_min, omega_ziel_vek[0])))
        self.omega[1] = np.sqrt(min(self.omega_quadrat_max, max(self.omega_quadrat_min, omega_ziel_vek[1])))
        self.omega[2] = np.sqrt(min(self.omega_quadrat_max, max(self.omega_quadrat_min, omega_ziel_vek[2])))
        self.omega[3] = np.sqrt(min(self.omega_quadrat_max, max(self.omega_quadrat_min, omega_ziel_vek[3])))

    def R(self):
        # siehe Gleichung 23
        return np.array([[cos(self.psi)*cos(self.theta), cos(self.psi)*sin(self.theta)*sin(self.phi) - sin(self.psi)*cos(self.phi), sin(self.psi)*sin(self.phi) + cos(self.psi)*sin(self.theta)*cos(self.phi)],
                         [sin(self.psi)*cos(self.theta), cos(self.psi)*cos(self.phi) + sin(self.psi)*sin(self.theta)*sin(self.phi), sin(self.psi)*sin(self.theta)*cos(self.phi) - cos(self.psi)*sin(self.phi)],
                         [-sin(self.theta),              cos(self.theta)*sin(self.phi),                                             cos(self.theta)*cos(self.phi)]])


    def welt_lineare_beschleunigungen(self):
        # siehe Gleichung 24
        R = self.R()
        g = np.array([0,0,self.g]).T
        c = -self.f_gesamt
        return g + np.matmul(R, np.array([0,0,c]).T) / self.m

    def koerper_achsen_beschleunigungen(self):
        # siehe Gleichung 27
        p_punkt = (self.m_x - self.r * self.q *(self.i_z - self.i_y))/self.i_x
        q_punkt = (self.m_y - self.r * self.p *(self.i_x - self.i_z))/self.i_y
        r_punkt = (self.m_z - self.q * self.p *(self.i_y - self.i_x))/self.i_z

        return np.array([p_punkt,q_punkt,r_punkt])


    def welt_achsen_beschleunigungen(self):
        # siehe Gleichung 28
        euler_rotationsmatrix = np.array([[1,         sin(self.phi) * tan(self.theta),  cos(self.phi) * tan(self.theta)],
                                          [0,         cos(self.phi),                   -sin(self.phi)],
                                          [0,         sin(self.phi) / cos(self.theta),  cos(self.phi) / cos(self.theta)]])

        return np.matmul(euler_rotationsmatrix, np.array([self.p,self.q,self.r]).T)

    def bringe_zustand_voran(self, dt):
        # siehe Gleichung 29
        wel_achs_bes = self.welt_achsen_beschleunigungen()
        koerp_achs_bes = self.koerper_achsen_beschleunigungen()
        welt_lin_bes = self.welt_lineare_beschleunigungen()

        X_punkt = np.array([self.X[6],
                            self.X[7],
                            self.X[8],
                            wel_achs_bes[0],
                            wel_achs_bes[1],
                            wel_achs_bes[2],
                            welt_lin_bes[0],
                            welt_lin_bes[1],
                            welt_lin_bes[2],
                            koerp_achs_bes[0],
                            koerp_achs_bes[1],
                            koerp_achs_bes[2]])

        self.X = self.X + X_punkt * dt      
        
class KaskadierendeRegelung3D():
    
    def __init__(self, z_k_p, z_k_d, x_k_p, x_k_d, y_k_p, y_k_d, k_p_roll, k_p_nick, k_p_gier, k_p_p, k_p_q, k_p_r, m, i_x, i_y, i_z):
        self.z_k_p = z_k_p
        self.z_k_d = z_k_d
        self.x_k_p = x_k_p
        self.x_k_d = x_k_d
        self.y_k_p = y_k_p
        self.y_k_d = y_k_d
        self.k_p_roll = k_p_roll
        self.k_p_nick = k_p_nick
        self.k_p_gier = k_p_gier
        self.k_p_p = k_p_p
        self.k_p_q = k_p_q
        self.k_p_r = k_p_r
        
        self.m = m
        self.i_x = i_x
        self.i_y = i_y
        self.i_z = i_z
        
        # Vermutlich haben die Ingenieure von der NASA deshalb eine Drohne auf den Mars geschickt, weil sie diese Zahl mal ändern wollten.
        # https://en.wikipedia.org/wiki/Ingenuity_(helicopter)
        self.g = 9.81
    
    
    def hoehen_regler(self, z_ziel, z_punkt_ziel, z_punkt_punkt_vorschub, z_aktuell, z_punkt_aktuell, rotationsmatrix):
        # siehe Gleichung 33
        z_punkt_punkt_ziel = self.z_k_p * (z_ziel - z_aktuell) + self.z_k_d * (z_punkt_ziel - z_punkt_aktuell) + z_punkt_punkt_vorschub
        
        # siehe Gleichung 34
        f_gesamt_ziel = -(self.m/rotationsmatrix[2,2])*(z_punkt_punkt_ziel - self.g)

        return f_gesamt_ziel


    def seiten_regler(self, x_ziel, x_punkt_ziel, x_punkt_punkt_vorschub, x_aktuell, x_punkt_aktuell, y_ziel, y_punkt_ziel, y_punkt_punkt_vorschub, y_aktuell, y_punkt_aktuell, f_gesamt_ziel):
        # siehe Gleichung 35
        x_punkt_punkt_ziel = self.x_k_p * (x_ziel - x_aktuell) + self.x_k_d * (x_punkt_ziel - x_punkt_aktuell) + x_punkt_punkt_vorschub
        y_punkt_punkt_ziel = self.y_k_p * (y_ziel - y_aktuell) + self.y_k_d * (y_punkt_ziel - y_punkt_aktuell) + y_punkt_punkt_vorschub
        
        # siehe Gleichung 36
        r_1_3_ziel = -(self.m*x_punkt_punkt_ziel)/f_gesamt_ziel
        r_2_3_ziel = -(self.m*y_punkt_punkt_ziel)/f_gesamt_ziel

        return r_1_3_ziel, r_2_3_ziel


    def roll_nick_regler(self, r_1_3_ziel, r_2_3_ziel, rotationsmatrix):            
        # siehe Gleichung 37
        r_1_3_ziel_punkt = self.k_p_roll * (r_1_3_ziel - rotationsmatrix[0,2])
        r_2_3_ziel_punkt = self.k_p_nick * (r_2_3_ziel - rotationsmatrix[1,2])

        # siehe Gleichung 38
        rotationsmatrix1 = np.array([[rotationsmatrix[1,0], -rotationsmatrix[0,0]], [rotationsmatrix[1,1], -rotationsmatrix[0,1]]]) / rotationsmatrix[2,2]
        pq_ziel = np.matmul(rotationsmatrix1, np.array([r_1_3_ziel_punkt, r_2_3_ziel_punkt]).T)
        p_ziel = pq_ziel[0]
        q_ziel = pq_ziel[1]

        return p_ziel, q_ziel


    def gier_regler(self, psi_ziel, psi_aktuell):
        # siehe Gleichung 39
        r_ziel = self.k_p_gier * (psi_ziel - psi_aktuell)

        return r_ziel


    def koerper_raten_regler(self, p_ziel, q_ziel, r_ziel, p_aktuell, q_aktuell, r_aktuell):
        # siehe Gleichung 40
        p_ziel_punkt = self.k_p_p * (p_ziel - p_aktuell)
        q_ziel_punkt = self.k_p_q * (q_ziel - q_aktuell)
        r_ziel_punkt = self.k_p_r * (r_ziel - r_aktuell)
        
        # siehe Gleichung 41
        m_x_ziel = self.i_x * p_ziel_punkt
        m_y_ziel = self.i_y * q_ziel_punkt
        m_z_ziel = self.i_z * r_ziel_punkt

        return m_x_ziel, m_y_ziel, m_z_ziel


# --- Simulationskonstanten ---
# siehe Kapitel "Ermittlung der Simulationsparameter"
k_f = 0.00141446535      
k_m = 0.0004215641
m = 1.0
L = 0.23 # siehe Abbildung "3D Drohne von oben" (Abbildungsverzeichnis)
i_x = 0.0121
i_y = 0.0119
i_z = 0.0223
omega_quadrat_min = 400
omega_quadrat_max = 4356

gesamte_simulationsdauer = 20.0
dt = 0.01
verhaeltnis_innere_aeussere_schleife = 10
# -----------------------------------------------------------------------------------

# --- Streckeninformationen ---
# siehe Kapitel "Dreidimensionale Simulation", Unterkapitel "Regelung", Paragraph "Strecke"
t = np.linspace(0.0, gesamte_simulationsdauer, int(gesamte_simulationsdauer/dt))

faktor_s_x = 0.8
faktor_s_y = 0.4
faktor_s_z = 0.4

s_x =  np.sin(faktor_s_x * t) 
s_x_punkt =  faktor_s_x * np.cos(faktor_s_x * t)
s_x_punkt_punkt = faktor_s_x**2 * -np.sin(faktor_s_x * t)

s_y =  np.cos(faktor_s_y * t)
s_y_punkt = faktor_s_y * -np.sin(faktor_s_y * t)
s_y_punkt_punkt = faktor_s_y**2 * -np.cos(faktor_s_y * t)

s_z = np.cos(faktor_s_z * t)
s_z_punkt = faktor_s_z * -np.sin(faktor_s_z * t)
s_z_punkt_punkt = faktor_s_z**2 * -np.cos(faktor_s_z * t)

s_psi = np.arctan2(s_y_punkt, s_x_punkt)
# -----------------------------------------------------------------------------------


# --- Simulation ---
drohne = Drohne3D(k_f, k_m, m, L, i_x, i_y, i_z, omega_quadrat_min, omega_quadrat_max)
drohne.X = np.array([s_x[0], s_y[0], s_z[0], 0.0, 0.0, s_psi[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
regelung = KaskadierendeRegelung3D(20.0, 4.25, 8.0, 4.0, 8.0, 4.0, 9.0, 9.0, 7.75, 22.0, 22.0, 2.5, m, i_x, i_y, i_z)

ax = plt.axes(projection='3d')
ax.plot(s_x, s_y, s_z, linestyle='-', marker='.', color='red')
punkte = []

# siehe Abbildung "Kaskadierende Regelarchitektur in drei Dimensionen" (Abbildungsverzeichnis)
for i in range(0, s_z.shape[0]):
    f_gesamt_ziel = regelung.hoehen_regler(s_z[i], s_z_punkt[i], s_z_punkt_punkt[i], drohne.X[2], drohne.X[8], drohne.R())
    r_1_3_ziel, r_2_3_ziel = regelung.seiten_regler(s_x[i], s_x_punkt[i], s_x_punkt_punkt[i], drohne.X[0], drohne.X[6], s_y[i], s_y_punkt[i], s_y_punkt_punkt[i], drohne.X[1], drohne.X[7], f_gesamt_ziel) 
    for _ in range(verhaeltnis_innere_aeussere_schleife):
        p_ziel, q_ziel = regelung.roll_nick_regler(r_1_3_ziel, r_2_3_ziel, drohne.R())
        r_ziel = regelung.gier_regler(s_psi[i], drohne.psi)
        m_x_ziel, m_y_ziel, m_z_ziel = regelung.koerper_raten_regler(p_ziel, q_ziel, r_ziel, drohne.X[9], drohne.X[10], drohne.X[11])
        drohne.setze_rotationsgeschwindigkeiten(f_gesamt_ziel, m_x_ziel, m_y_ziel, m_z_ziel)
        drohne.bringe_zustand_voran(dt/verhaeltnis_innere_aeussere_schleife)


    rotationsmatrix = drohne.R()
    # Drohnenkörpermittelpunkt
    punkte.append(ax.scatter3D(drohne.x, drohne.y, drohne.z, color='blue'))
    f = [drohne.f_1, drohne.f_2, drohne.f_3, drohne.f_4]
    zaehler = 0
    
    for j in range(-1, 2, 2):
        for k in range(-1, 2, 2):
            # Propeller
            prop_pos = np.matmul(rotationsmatrix, np.array([drohne.l*j, drohne.l*k, 0]).T)
            punkte.append(ax.scatter3D(drohne.x+prop_pos[0], drohne.y+prop_pos[1], drohne.z+prop_pos[2], color='green'))
            
            # Kraftvektoren
            kraft_pos = np.matmul(rotationsmatrix, np.array([drohne.l*j, drohne.l*k, 0-f[zaehler]]).T)
            punkte.append(ax.scatter3D(drohne.x+kraft_pos[0], drohne.y+kraft_pos[1], drohne.z+kraft_pos[2], color='red'))
            
            zaehler+=1

    plt.pause(0.01)
    
    for p in punkte:
        p.remove()

    punkte.clear()

# -----------------------------------------------------------------------------------