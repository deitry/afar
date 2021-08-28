! = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! Вспомогательный модуль расчёта коэффициентов и др.
!
! Ворнычев Д.С., асп. МГТУ им.Баумана, каф.Э6
! Дата последнего изменения: 2014.05.14
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

MODULE FUNC

USE PROG_DATA

IMPLICIT NONE


CONTAINS
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Вычисление коэффициента температуропроводности для точки с 
! данными координатами
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL FUNCTION ATEMP(X, Y, Z) RESULT(res)

REAL :: X, Y, Z

	IF (Z <= delta1) THEN
		res = AT1
		RETURN
	END IF

	res = AT2

END FUNCTION ATEMP

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! расчёт коэффициента теплопроводности
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL FUNCTION LAMBDA(X,Y,Z) RESULT(res)
	REAL	::	X,Y,Z
	
	IF (Z <= delta1) THEN
		res = LT1
		RETURN
	END IF

	res = LT2

END FUNCTION LAMBDA

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Вычисление внутренних источников
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
REAL FUNCTION QVN(X, Y, Z) RESULT(res)
	! TODO : добавить T в данной точке в качестве параметра? Можно будет
	! учитывать нелинейности. Но пока их нет - так тоже ничего

	REAL 			:: 	X, Y, Z
	INTEGER 		:: 	m
	TYPE(POWEL) 	:: 	p

	! если отсутствуют тепловыделяющие элементы
	IF (pCnt == 0) THEN
		res = 0
		RETURN
	END IF

	! проходимся по списку "мощностей" и смотрим, не попали ли. Если попали -
	! возвращаем значение мощности
	DO m = 0, pCnt-1

		p = powers(m)
		
		! если наши координаты попадают в пределы x0<x<(x0+lx) и y0<y<(y0+ly)
		! считаем, что тепловыделяющие элементы присутствуют только на границе
		IF (		(X >= p%x0).AND.(X <= (p%x0 + p%lx))	&
	&		.AND.	(Y >= p%y0).AND.(Y <= (p%y0 + p%ly))	&
	&		.AND. 	(Z == 0)	) 	THEN
		
			! то возвращаем значение мощности для данного элемента
			res = p%Pow
			RETURN
		END IF
	END DO
		
	res = 0

END FUNCTION QVN

REAL FUNCTION TZM(k) RESULT(res)
    INTEGER :: i,j,k
    
    res = 0
    
    DO i = 0, Mi
        DO j = 0, Mj
            res = res + T(i,j,k,0)/((Mi+1)*(Mj+1))
        END DO
    END DO

END FUNCTION TZM

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ВЫВОД ВСЕГО СЛОЯ В ФАЙЛ
SUBROUTINE PRINT_FIELD(fname)
	
	CHARACTER*512 :: fname, fmt
	INTEGER :: unum = 33, jy, jk

	OPEN(UNIT=unum,FILE=TRIM(fname), STATUS='REPLACE')

    ! средняя температура по слоям
    DO jk = 0, Mk
        WRITE(unum, '(f6.3, a,f10.4)') jk*hz, '	', TZM(jk) 
    END DO

	! цикл по Z
	DO jk = 0, Mk

		WRITE(unum, *) "Слой Z = ", jk*hz

		DO jy = 0, Mj
			WRITE(fmt,'( "(", I0,"(F0.0,",a,"))" )') Mi+1, "'	'"
			!WRITE(*,*) fmt
			WRITE(unum,fmt) T(:,jy,jk,0)
		END DO

		WRITE(unum, *)
		WRITE(unum, *)
		
	END DO

	CLOSE(unum)

END SUBROUTINE PRINT_FIELD


END MODULE FUNC

