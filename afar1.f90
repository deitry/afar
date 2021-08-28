! = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
! Основной модуль численного расчёта уравнения теплопроводности
! для двуслойной мат. модели ППМ АФАР.
!
! Текущий вариант решения: нестационарная двумерная задача, 
! неявная схема с разделением по направлениям.
!
! Ворнычев Д.С., асп. МГТУ им.Баумана, каф.Э6
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

PROGRAM AFAR1

USE FUNC
USE PROG_DATA
USE SCHEME3

IMPLICIT NONE

CALL FILE_NAMES()

! открытие канала для ведения лога
OPEN (UNIT=output , FILE=outFile , STATUS='REPLACE')

WRITE(*, *) "Initialization... "
WRITE(output, *) "Инициализация... "
	
! инициализация данных
CALL DATA_INIT()

!WRITE(*, *) "Успех!"
!WRITE(output, *) "Успех!"
	

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! цикл расчёта по времени
WRITE(*,*) "Calculation start..."
iTime: DO n = 0, Mt

	CALL NEXT_STEP()
	
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	! TODO : промежуточный вывод данных
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	IF (mod(n,INT(tprint/tau)) == 0) THEN
		WRITE(*,*) "Time = ", n*tau
		WRITE(out2File, '(a,F9.4, a )') timeFold//delimiter, n*tau, ext	!REAL(INT(n*tau*1000))/1000
		CALL PRINT_FIELD(out2File)
	END IF

END DO iTime

!WRITE(*, *) "Деаллокация массивов, закрытие каналов"
!WRITE(output, *) "Деаллокация массивов, закрытие каналов"
	

! закрытие каналов, высвобождение памяти
CALL DATA_FINALIZE()

WRITE(*, *) "End"
WRITE(output, *) "End"
	
CLOSE(output)

END PROGRAM AFAR1

