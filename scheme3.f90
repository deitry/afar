! Cюда будет вынесена схема решения дифура в виде функции 
! для перехода на следующий шаг по времени

! Решение трёхмерного уравнения теплопроводности

MODULE SCHEME3

USE PROG_DATA
USE FUNC

IMPLICIT NONE

! вспомогательные переменные
REAL 		:: 	Ai, Bi, Ci, Fi

PRIVATE 	:: 	Ai, Bi, Ci, Fi

CONTAINS

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! функции для расчёта вспомогательных коэффициентов 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! - - - - - - - - - - - - - - - - - - - - - - - - -
! для первого этапа интегрирования
REAL FUNCTION Ai1(i,j,k)	! при i,j-1,k
	INTEGER :: i,j,k
	Ai1 = tau/2 * ATEMP(i*hx, j*hy, k*hz) / hx2
END FUNCTION Ai1

REAL FUNCTION Bi1(i,j,k)	! при i,j+1,k
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

! - - - - - - - - - - - - - - - - - - - - - - - - -
! для второго этапа интегрирования
REAL FUNCTION Ai2(i,j,k)	! при i,j-1,k
	INTEGER :: i,j,k
	Ai2 = tau/2 * ATEMP(i*hx, j*hy, k*hz) / hy2
END FUNCTION Ai2

REAL FUNCTION Bi2(i,j,k)	! при i,j+1,k
	INTEGER :: i,j,k
	Bi2 = tau/2 * ATEMP(i*hx, j*hy, k*hz) / hy2
END FUNCTION Bi2

REAL FUNCTION Ci2(i,j,k)	! при i,j,k
	INTEGER :: i,j,k
	Ci2 = - tau * ATEMP(i*hx, j*hy, k*hz) / hy2 - 1
END FUNCTION Ci2

REAL FUNCTION Fi2(i,j,k)	! правая часть
	INTEGER :: i,j,k
	Fi2 = - T(i,j,k,1)
END FUNCTION Fi2

! - - - - - - - - - - - - - - - - - - - - - - - - -
! для третьего этапа
REAL FUNCTION Ai3(i,j,k)	! при i,j-1,k
	INTEGER :: i,j,k
	Ai3 = tau/2 * ATEMP(i*hx, j*hy, k*hz) / hz2
END FUNCTION Ai3

REAL FUNCTION Bi3(i,j,k)	! при i,j+1,k
	INTEGER :: i,j,k
	Bi3 = tau/2 * ATEMP(i*hx, j*hy, k*hz) / hz2
END FUNCTION Bi3

REAL FUNCTION Ci3(i,j,k)	! при i,j,k
	INTEGER :: i,j,k
	Ci3 = - tau * ATEMP(i*hx, j*hy, k*hz) / hz2 - 1
END FUNCTION Ci3

REAL FUNCTION Fi3(i,j,k)	! правая часть
	INTEGER :: i,j,k
	Fi3 = - T(i,j,k,2)
END FUNCTION Fi3

! - - - - - - - - - - - - - - - - - - - - - - - - -
! четвёртый этап интегрирования - "корректор"
REAL FUNCTION Prav4(i,j,k)
	INTEGER ::	i,j,k
	REAL 	::	L1, L2, L3
	REAL	::	As, Qs 

	IF (i == 0) THEN
		L1 = (T(i+1,j,k,3) - T(i,j,k,3)) / hx2
	ELSE IF (i == Mi) THEN
		L1 = (T(i-1,j,k,3) - T(i,j,k,3)) / hx2
	ELSE
		L1 = (T(i-1,j,k,3) - 2*T(i,j,k,3) + T(i+1,j,k,3)) / hx2
	END IF
	
	IF (j == 0) THEN
		L2 = (T(i,j+1,k,3) - T(i,j,k,2)) / hy2
	ELSE IF (j == Mj) THEN
		L2 = (T(i,j-1,k,3) - T(i,j,k,2)) / hy2
	ELSE
		L2 = (T(i,j-1,k,3) - 2*T(i,j,k,2) + T(i,j+1,k,2)) / hy2
	END IF
	
	IF (k == 0) THEN
		L3 = (T(i,j,k+1,3) - T(i,j,k,3)) / hz2
	ELSE IF (k == Mk) THEN
		L3 = (T(i,j,k-1,3) - T(i,j,k,3)) / hz2
	ELSE
		L3 = (T(i,j,k-1,3) - 2*T(i,j,k,3) + T(i,j,k+1,3)) / hz2
	END IF
	
	!Prav4 = T(i,j,k,3)
	As = ATEMP	(i*hx, j*hy, k*hz)
	Qs = QVN	(i*hx, j*hy, k*hz)
	Prav4 = T(i,j,k,0) + tau*As * (L1+L2+L3) + tau*Qs
	END FUNCTION Prav4

! =	= = = = = = = = = = = = = = = = = = = = = = = = = = = = 
! основной расчёт
! =	= = = = = = = = = = = = = = = = = = = = = = = = = = = = 

! переход на следующий шаг расчёта по времени
SUBROUTINE NEXT_STEP()
	
	INTEGER :: i, j, k, kt
	
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	! этап 1: первый промежуточный слой; прогонка по X
	E1: DO k = 0, Mk
		DO j = 0, Mj
		
		! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		! - вычисление прогоночных коэффициентов
		!alphaX(1) = - Bix(0,j)/Cix(0,j)
		!betaX(1) = Fix(0,j)/Cix(0,j)
		! при нулевом потоке через границу
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
		!T(Mi,j,1) = (Fix(Mi,j) - Aix(Mi,j)*betaX(Mi))/(Cix(Mi,j) + Aix(Mi,j)*alphaX(Mi))
		T(Mi,j,k,1) = beta1(Mi) / (1 - alpha1(Mi))	! при нулевом потоке через границу
		
		E1S2: DO i = Mi-1, 0, -1
			T(i,j,k,1) = alpha1(i+1)*T(i+1,j,k,1) + beta1(i+1)
		END DO E1S2
	
	END DO		! j
	END DO E1 	! k
	
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	! этап 2: прогонка по Y
	E2: DO i = 0, Mi
		DO k = 0, Mk 
		! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		! - вычисление прогоночных коэффициентов
		!alphaY(1) = - Biy(i,0)/Ciy(i,0)
		!betaY(1) = Fiy(i,0)/Ciy(i,0)
		alpha2(1) = 1
		beta2(1) = 0
		
		! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		E2S1: DO j = 1, Mj-1
			! многократно встречающиеся коэффициенты
			Ai = Ai2(i,j,k) 		!к-т при Т(i-1,j)
			Ci = Ci2(i,j,k)			!к-т при Т(i,j)
			
			alpha2(j+1) = - Bi2(i,j,k) / (Ai*alpha2(j)+Ci)
			beta2(j+1) = (Fi2(i,j,k) - Ai*beta2(j)) / (Ai*alpha2(j)+Ci)

		END DO E2S1
	
		! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		! - вычисление окончательных значений
		T(i,Mj,k,2) = beta2(Mj) / (1-alpha2(Mj))
		
		E2S2: DO j = Mj-1, 0, -1
			T(i,j,k,2) = alpha2(j+1)*T(i,j+1,k,2) + beta2(j+1)
		END DO E2S2
		
	END DO
	END DO E2
	
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	! этап 3: третий промежуточный слой; решение по Z. Решается прогонкой
	E3: DO i = 0, Mi
		DO j = 0, Mj
		
		! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		! - вычисление прогоночных коэффициентов
		As = LAMBDA(i*hx,j*hy,Mk*hz)
		
		alpha3(1) = As/(As + hz*KTA1)
		beta3(1)  = hz*KTA1/(As + hx*KTA1) * Te1
		
		E3S1: DO k = 1, Mk

			Ai = Ai3(i,j,k)
			Ci = Ci3(i,j,k)
			
			alpha3(k+1) = - Bi3(i,j,k) / (Ai*alpha3(k)+Ci)
			beta3(k+1) = (Fi3(i,j,k) - Ai*beta3(k)) / (Ai*alpha3(k)+Ci)
		
		END DO E3S1
		
		! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		! - вычисление промежуточных значений
		!T(i,j,Mk,3) = (Fi3(i,j,Mk) - Ai3(i,j,Mk)*beta3(Mk))/(Ci3(i,j,Mk) + Ai3(i,j,Mk)*alpha3(Mk))
		!T(i,j,Mk,3) = beta3(Mk) / (1 - alpha3(Mk))	! при нулевом потоке через границу
		!T(i,j,Mk,3)	= (As*beta3(Mk) - qk2) / (As*(1-alpha3(Mk)))
		
		As = LAMBDA(i*hx,j*hy,Mk*hz)
		T(i,j,Mk,3) = (As*beta3(Mk) + hz*KTA2*Te2) / (hz*KTA2 + As*(1 - alpha3(Mk)))
		
		
		E3S2: DO k = Mk-1, 0, -1
			T(i,j,k,3) = alpha3(k+1)*T(i,j,k+1,3) + beta3(k+1)
		END DO E3S2
	
	END DO		! j
	END DO E3 	! 

	! четвёртый этап интегрирования - "корректор"
	DO i = 0, Mi
		DO j = 0, Mj
			DO k = 0, Mk
				T(i,j,k,4) = Prav4(i,j,k)
			END DO
		END DO
	END DO

	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	! переход на следующий уровень
	CALL MOVE_STEP()
	
END SUBROUTINE NEXT_STEP

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! переход на следующий временной слой
SUBROUTINE MOVE_STEP()
	DO i = 0, Mi
		DO j = 0, Mj
			DO k = 0, Mk
				T(i,j,k,0) = T(i,j,k,layCnt)		! обновляем текущее
			END DO
		END DO
	END DO
END SUBROUTINE MOVE_STEP

END MODULE SCHEME3

