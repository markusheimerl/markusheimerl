/*
* Markus Heimerl
* 2021, OTH Regensburg
* markus (at) markusheimerl (dot) com
*/

#include <math.h>

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef EINGRENZEN
#define EINGRENZEN(a,untere_schranke,obere_schranke) MAX((untere_schranke),MIN((a),(obere_schranke)))
#endif

#define F_PI ((float)(M_PI))

			
// TODO: Stelle sicher, dass die Rotationsrichtungen der Motoren wie folgt verbaut sind:
// Vorne links und hinten rechts im Uhrzeigersinn
// Vorne rechts und hinten links gegen den Uhrzeigersinn

class KaskadierendeRegelung3D{
	
	public:
		KaskadierendeRegelung3D(){
			g = 9.81f;
			integrierterHoehenfehler = 0.0; // Diese Variable läuft irgendwann über. Testflüge sollten kurz gehalten werden.
			dt = 0.0018; // Der Regler muss also exakt alle 0.0018 Sekunden (oder 1.8 ms) aufgerufen werden.
			
			m = 1.0;
			L = 0.23;
			l = (L * sqrt(2)) / 2;
			i_x = 0.0121;
			i_y = 0.0119;
			i_z = 0.0223;
			
			k_f = 0.00141446535;
			k_m = 0.0004215641;
			
			omega_quadrat_min = 400;
			omega_quadrat_max = 4356;
			
			maxSinkrate = 2.0;
			maxSteigrate = 5.0;
			maxLateralBeschl = 5.0;
			maxLageWinkel = 0.7;
			minSchub = 0.56898;
			maxSchub = 6.0822;
			
			z_k_p = 20.0;
			z_k_d = 4.25;
			z_k_i = 5.0;
			x_k_p = 8.0;
			x_k_d = 4.0;
			y_k_p = 8.0;
			y_k_d = 4.0;
			k_p_roll = 9.0;
			k_p_nick = 9.0;
			k_p_gier = 7.5;
			k_p_p = 22.0;
			k_p_q = 22.0;
			k_p_r = 2.5;
		}
	
		// Hauptfunktion; Diese muss von außen aufgerufen werden.
		void regle(float x_aktuell, float y_aktuell, float z_aktuell, float phi_aktuell, float theta_aktuell, float psi_aktuell, float p_aktuell, float q_aktuell, float r_aktuell, 
					float x_punkt_aktuell, float y_punkt_aktuell, float z_punkt_aktuell, float x_punkt_punkt_vorschub, float y_punkt_punkt_vorschub, float z_punkt_punkt_vorschub,
						float x_ziel, float y_ziel, float z_ziel, float psi_ziel, float x_punkt_ziel, float y_punkt_ziel, float z_punkt_ziel, 
							float* pwm_1, float* pwm_2, float* pwm_3, float* pwm_4){
			
			
			float rotationsmatrix[3][3];
			
			
			hole_rotationsmatrix(phi_aktuell, theta_aktuell, psi_aktuell, rotationsmatrix);
			
			
			float f_gesamt_ziel = hoehen_regler(z_ziel, z_punkt_ziel, z_punkt_punkt_vorschub, z_aktuell, z_punkt_aktuell, rotationsmatrix);
			// Spare ein bisschen Schub für die Winkelkontrolle auf
			float schubSpielraum = 0.1f * (maxSchub - minSchub);
			f_gesamt_ziel = EINGRENZEN(f_gesamt_ziel, (minSchub + schubSpielraum) * 4.f, (maxSchub - schubSpielraum) * 4.f);
			
			
			float r_1_3_ziel, r_2_3_ziel;
			seiten_regler(x_ziel, x_punkt_ziel, x_punkt_punkt_vorschub, x_aktuell, x_punkt_aktuell, y_ziel, y_punkt_ziel, y_punkt_punkt_vorschub, y_aktuell, y_punkt_aktuell, f_gesamt_ziel, &r_1_3_ziel, &r_2_3_ziel);
			
			
			float p_ziel, q_ziel;
			roll_nick_regler(r_1_3_ziel, r_2_3_ziel, rotationsmatrix, &p_ziel, &q_ziel);
			
			
			float r_ziel = gier_regler(psi_ziel, psi_aktuell);
			
			
			float m_x_ziel, m_y_ziel, m_z_ziel;
			koerper_raten_regler(p_ziel, q_ziel, r_ziel, p_aktuell, q_aktuell, r_aktuell, &m_x_ziel, &m_y_ziel, &m_z_ziel);
			
			
			float omega_1, omega_2, omega_3, omega_4;
			hole_rotationsgeschwindigkeiten(f_gesamt_ziel, m_x_ziel, m_y_ziel, m_z_ziel, &omega_1, &omega_2, &omega_3, &omega_4);
			
			
			hole_pwm_signal(omega_1, omega_2, omega_3, omega_4, pwm_1, pwm_2, pwm_3, pwm_4);
		}
	
	
	private:
		float hoehen_regler(float z_ziel, float z_punkt_ziel, float z_punkt_punkt_vorschub, float z_aktuell, float z_punkt_aktuell, float rotationsmatrix[3][3]){
			
			float z_postition_fehler = z_ziel - z_aktuell;
			integrierterHoehenfehler += z_postition_fehler * dt;
			
			// siehe Gleichung 33
			float z_punkt_punkt_ziel = z_k_p * z_postition_fehler + z_k_d * (z_punkt_ziel - z_punkt_aktuell) + z_k_i * integrierterHoehenfehler + z_punkt_punkt_vorschub;
			
			// limitiere Höhenbeschleunigung
			z_punkt_punkt_ziel = EINGRENZEN(z_punkt_punkt_ziel, -maxSinkrate / dt, maxSteigrate / dt);
			
			// siehe Gleichung 34
			float f_gesamt_ziel = -(m/rotationsmatrix[2][2])*(z_punkt_punkt_ziel - g);
			
			// Lass keine gewollte Kraft kleiner 0 zu; Die Propeller können keine Kraft nach unten wirken und der Nenner darf im Seitenregler nicht null werden.
			if(f_gesamt_ziel <= 0.0) f_gesamt_ziel = 0.0001f;

			return f_gesamt_ziel;
		}

		void seiten_regler(float x_ziel, float x_punkt_ziel, float x_punkt_punkt_vorschub, float x_aktuell, float x_punkt_aktuell, float y_ziel, float y_punkt_ziel, float y_punkt_punkt_vorschub, 
								float y_aktuell, float y_punkt_aktuell, float f_gesamt_ziel, float* r_1_3_ziel, float* r_2_3_ziel){
			
			// siehe Gleichung 35
			float x_punkt_punkt_ziel = x_k_p * (x_ziel - x_aktuell) + x_k_d * (x_punkt_ziel - x_punkt_aktuell) + x_punkt_punkt_vorschub;
			float y_punkt_punkt_ziel = y_k_p * (y_ziel - y_aktuell) + y_k_d * (y_punkt_ziel - y_punkt_aktuell) + y_punkt_punkt_vorschub;
			
			// limitiere Lateralbeschleunigung
			x_punkt_punkt_ziel = EINGRENZEN(x_punkt_punkt_ziel, -maxLateralBeschl / dt, maxLateralBeschl / dt);
			y_punkt_punkt_ziel = EINGRENZEN(y_punkt_punkt_ziel, -maxLateralBeschl / dt, maxLateralBeschl / dt);
			
			// siehe Gleichung 36
			*r_1_3_ziel = -(m*x_punkt_punkt_ziel)/f_gesamt_ziel;
			*r_2_3_ziel = -(m*y_punkt_punkt_ziel)/f_gesamt_ziel;
		}

		void roll_nick_regler(float r_1_3_ziel, float r_2_3_ziel, float rotationsmatrix[3][3], float* p_ziel, float* q_ziel){
		
			// limitiere Lagewinkel
			r_1_3_ziel = EINGRENZEN(r_1_3_ziel, -maxLageWinkel, maxLageWinkel);
			r_2_3_ziel = EINGRENZEN(r_2_3_ziel, -maxLageWinkel, maxLageWinkel);
		
			// siehe Gleichung 37
			float r_1_3_ziel_punkt = k_p_roll * (r_1_3_ziel - rotationsmatrix[0][2]);
			float r_2_3_ziel_punkt = k_p_nick * (r_2_3_ziel - rotationsmatrix[1][2]);

			// siehe Gleichung 38
			*p_ziel = rotationsmatrix[1][0] * r_1_3_ziel_punkt + (-rotationsmatrix[0][0]) * r_2_3_ziel_punkt;
			*q_ziel = rotationsmatrix[1][1] * r_1_3_ziel_punkt + (-rotationsmatrix[0][1]) * r_2_3_ziel_punkt;
		}
		
		float gier_regler(float psi_ziel, float psi_aktuell){
			
			// Korriegere Psi-Ziel
			float psi_ziel_korrigiert = 0;
			if(psi_ziel > 0){
				psi_ziel_korrigiert = fmodf(psi_ziel, 2 * F_PI);
			}else{
				psi_ziel_korrigiert = -fmodf(psi_ziel, 2 * F_PI);
			}
			
			// siehe Gleichung 39
			float r_ziel = k_p_gier * (psi_ziel_korrigiert - psi_aktuell);

			// Limitiere r_ziel (schließe Kreis)
			if(r_ziel > F_PI){
				r_ziel -= 2 * F_PI;
			}
			
			if(r_ziel < -F_PI){
				r_ziel += 2 * F_PI;
			}

			return r_ziel;
		}

		void koerper_raten_regler(float p_ziel, float q_ziel, float r_ziel, float p_aktuell, float q_aktuell, float r_aktuell, float* m_x_ziel, float* m_y_ziel, float* m_z_ziel){
			// siehe Gleichung 40
			float p_ziel_punkt = k_p_p * (p_ziel - p_aktuell);
			float q_ziel_punkt = k_p_q * (q_ziel - q_aktuell);
			float r_ziel_punkt = k_p_r * (r_ziel - r_aktuell);
			
			// siehe Gleichung 41
			*m_x_ziel = i_x * p_ziel_punkt;
			*m_y_ziel = i_y * q_ziel_punkt;
			*m_z_ziel = i_z * r_ziel_punkt;
		}
		
		void hole_rotationsmatrix(float phi, float theta, float psi, float rotationsmatrix[3][3]){
			// siehe Gleichung 23
			rotationsmatrix[0][0] = cos(psi)*cos(theta);
			rotationsmatrix[0][1] = cos(psi)*sin(theta)*sin(phi) - sin(psi)*cos(phi);
			rotationsmatrix[0][2] = sin(psi)*sin(phi) + cos(psi)*sin(theta)*cos(phi);
			
			rotationsmatrix[1][0] = sin(psi)*cos(theta);
			rotationsmatrix[1][1] = cos(psi)*cos(phi) + sin(psi)*sin(theta)*sin(phi);
			rotationsmatrix[1][2] = sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi);
			
			rotationsmatrix[2][0] = -sin(theta);
			rotationsmatrix[2][1] = cos(theta)*sin(phi);
			rotationsmatrix[2][2] = cos(theta)*cos(phi);
		}
		
		void hole_rotationsgeschwindigkeiten(float f_gesamt_ziel, float m_x_ziel, float m_y_ziel, float m_z_ziel, float* omega_1, float* omega_2, float* omega_3, float* omega_4){
			// siehe Gleichung 44
			float ziel_vek[4] = {(f_gesamt_ziel/k_f), (m_x_ziel/(l*k_f)), (m_y_ziel/(l*k_f)), (m_z_ziel/(k_m))};
		
			float o1 = ziel_vek[0] * 0.25 + ziel_vek[1] * 0.25 + ziel_vek[2] * 0.25 + ziel_vek[3] * 0.25;
			*omega_1 = sqrtf(EINGRENZEN(o1, omega_quadrat_min, omega_quadrat_max));

			float o2 = ziel_vek[0] * 0.25 + ziel_vek[1] * -0.25 + ziel_vek[2] * 0.25 + ziel_vek[3] * -0.25;
			*omega_2 = sqrtf(EINGRENZEN(o2, omega_quadrat_min, omega_quadrat_max));
			
			float o3 = ziel_vek[0] * 0.25 + ziel_vek[1] * 0.25 + ziel_vek[2] * -0.25 + ziel_vek[3] * -0.25;
			*omega_3 = sqrtf(EINGRENZEN(o3, omega_quadrat_min, omega_quadrat_max));
			
			float o4 = ziel_vek[0] * 0.25 + ziel_vek[1] * -0.25 + ziel_vek[2] * -0.25 + ziel_vek[3] * 0.25;
			*omega_4 = sqrtf(EINGRENZEN(o4, omega_quadrat_min, omega_quadrat_max));
		}
		
		void hole_pwm_signal(float omega_1, float omega_2, float omega_3, float omega_4, float* pwm_1, float* pwm_2, float* pwm_3, float* pwm_4){
			// siehe 06_messurements -> fit_k_pwm.py und "rotorgeschwindigkeit v pwm.txt"
			*pwm_1 = 0.15277 * (omega_1*omega_1) + 3.5900 * omega_1 + 1080.91881;
			*pwm_2 = 0.15277 * (omega_2*omega_2) + 3.5900 * omega_2 + 1080.91881;
			*pwm_3 = 0.15277 * (omega_3*omega_3) + 3.5900 * omega_3 + 1080.91881;
			*pwm_4 = 0.15277 * (omega_4*omega_4) + 3.5900 * omega_4 + 1080.91881;
		}

		float integrierterHoehenfehler;
		float z_k_p, z_k_d, z_k_i, x_k_p, x_k_d, y_k_p, y_k_d, k_p_roll, k_p_nick, k_p_gier, k_p_p, k_p_q, k_p_r;
		float m, L, l;
		float k_f, k_m;
		float g;
		float dt;
		float i_x, i_y, i_z;
		float maxSinkrate, maxSteigrate;
		float maxLateralBeschl;
		float maxLageWinkel;
		float maxSchub, minSchub;
		float omega_quadrat_min, omega_quadrat_max;
};
