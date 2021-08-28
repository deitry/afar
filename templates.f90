! файл с различного рода шаблонами для дальнейшего использования в программах ВТГ

! ===================================================================================
! Алгоритм прогонки
! ===================================================================================
! вспомогательные коэффициенты
REAL FUNCTION Ai1(i,j,k)	! при i-1,j,k
	INTEGER :: i,j,k
	Ai1 = tau/2 * ATEMP(i*hx, j*hy, k*hz) / hx2
END FUNCTION Ai1

REAL FUNCTION Bi1(i,j,k)	! при i+1,j,k
	INTEGER :: i,j,k
	Bi1 = tau/2 * ATEMP(i*hx, j*hy, k*hz) / hx2
END FUNCTION Bi1

REAL FUNCTION Ci1(i,j,k)	! при i,j,k
	INTEGER :: i,j,k
	Ci1 = - tau * ATEMP(i*hx, j*hy, k*hz) / hx2 - 1
END FUNCTION Ci1

REAL FUNCTION Fi1(i,j,k)	! правая часть
	INTEGER :: i,j,k
	Fi1 = - T(i,j,k,0) - QVN(i*hx,j*hy,k*hz) * tau/2
END FUNCTION Fi1
! ===================================================================================
! основной алгоритм
E1: DO k = 0, Mk
		DO j = 0, Mj
		
		! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		! - вычисление прогоночных коэффициентов
		alpha1(1) = 1
		beta1(1) = 0
		
		E1S1: DO i = 1, Mi-1

			Ai = Ai1(i,j,k)
			Ci = Ci1(i,j,k)
			
			alpha1(i+1) = - Bi1(i,j,k) / (Ai*alpha1(i)+Ci)
			beta1(i+1) = (Fi1(i,j,k) - Ai*beta1(i))/(Ai*alpha1(i)+Ci)
		
		END DO E1S1
		
		! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		! - вычисление промежуточных значений
		T(Mi,j,k,1) = beta1(Mi) / (1 - alpha1(Mi))
		
		E1S2: DO i = Mi-1, 0, -1
			T(i,j,k,1) = alpha1(i+1)*T(i+1,j,k,1) + beta1(i+1)
		END DO E1S2
	
	END DO		! j
	END DO E1 	! k
	
	
! ===================================================================================
! Вычисление крайних значений в прогонке при различных краевых условиях
! ===================================================================================

! = = = = = = = = = = = = = = = = = = 
! ГУ 1го рода (T)
alpha1(1) = - Bi1(0,j,k)/Cix(0,j,k)
beta1(1) = Fi1(0,j,k)/Cix(0,j,k)

T(Mi,j,1) = (Fix(Mi,j) - Aix(Mi,j)*betaX(Mi))/(Cix(Mi,j) + Aix(Mi,j)*alphaX(Mi))

! = = = = = = = = = = = = = = = = = = 
! ГУ 2го рода (q)
! - первый порядок аппроксимации
alpha1(1) = 1

! - - - - 
! при q = 0
beta1(1) = 0			
T(Mi,j,k,1) = beta3(Mk) / (1 - alpha3(Mk))	! при q = 0

! - - - -
! при q != 0
beta1(1) = hx*qk1/LAMBDA	
T(Mi,j,k,1)	= (LAMBDA*beta1(Mi) - qk2) / (LAMBDA*(1-alpha1(Mi))) ! 

! - - - - - - - - - - 
! - второй порядок аппроксимации

! = = = = = = = = = = = = = = = = = =
! ГУ 3го рода (T+q)
! - первый порядок аппроксимации
alpha1(0) = LAMBDA/(LAMBDA + hx*KTA1)	! KTA - к-т теплоотдачи: KTA*(T - Te)
beta1(0)  = hx*KTA1/(LAMBDA + hx*KTA1) * Te1

T(Mi,j,k,1) = (LAMBDA*beta(Mi) + hx*KTA2*Te2) / (hx*KTA2 + LAMBDA*(1-alpha1(Mi)))
! - второй порядок аппроксимации
