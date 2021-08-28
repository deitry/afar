! Модуль, описывающий основные переменные, а также хранящий процедуру
! по чтению данных из файла
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MODULE PROG_DATA

IMPLICIT NONE

INTEGER 			:: 		OS_KEY=1					! 0 - linux, 1 - windows		

REAL, ALLOCATABLE, DIMENSION (:,:,:,:) :: T 			! искомое поле температур
														! 1,2,3 индекс - пространственные координаты
														! 4 индекс - шаг по времени, в т.ч. предыдущий

INTEGER 			::		i, j, k, 				&	! пространственные индексы
	&						n							! временной индекс	

INTEGER 			::		Mi=20, Mj=20, Mk=5			! количество точек - вычисляется
INTEGER 			:: 		Mt	! количество точек по времени - вычисляется

REAL				::		Tmax = 5, 				&	! время расчёта
	&						T0 = 0, 				& 	! начальное значение поля температур
	&						AT1, AT2,				& 	! коэффициенты температуропроводности "верхней" и "нижней" пластины
	&						hx = 0.1, 				& 	! шаг расчёта вдоль оси x
	&						hy = 0.1, 				& 	! шаг расчёта вдоль оси y
	&						hz = 0.1, 				& 	! шаг расчёта вдоль оси z
	&						tau = 0.0001, 			& 	! шаг расчёта по времени
	&						tprint = 0.5,			&	! шаг вывода на печать
	&						etaX=0.5, etaY=0.5			! к-ты "смешивания" явной и неявной схемы
	
REAL 				:: 		hx2, hy2, hz2				! hx2=hx**2

INTEGER, PARAMETER	::		layCnt = 4, 			&	! количество слоёв по времени 
&							prevCnt = 0					! сколько значений по времени надо сохранять минус один

INTEGER 			:: 		output = 10					! канал для вывода данных

REAL, ALLOCATABLE, DIMENSION (:) :: alpha1, beta1, &	! прогоночные коэффициенты
&									alpha2, beta2, &
&									alpha3, beta3 

REAL 				:: 		lx=1., ly=1., lz=1.			! габариты всей области
REAL 				:: 		delta1						! толщина "верхней" пластины. Толщина нижней = lz - delta1

! тип, определяющий "мощность" - POWer ELement
TYPE 				:: 		POWEL
	REAL 			:: 		x0, y0, lx, ly				! положение элемента и его габариты
	REAL			::		Pow							! выделяемая мощность
END TYPE POWEL

TYPE(POWEL), ALLOCATABLE, DIMENSION(:) :: powers		! массив, описывающий геометрию
INTEGER 			:: 		pCnt						! количество элементов

CHARACTER*512 		:: 		inFile, geomFile, outFile, out2File
CHARACTER 			:: 		delimiter, TAB = CHAR(9)
CHARACTER*4 		::		timeFold
CHARACTER*5 		:: 		ext
	
REAL				::		qk1 = 1000, qk2 = 1000	! теплоотвод от поверхностей
REAL				::		As, Qs	! для вспомогательных вычислений
REAL				::		LT1, LT2
REAL				::		KTA1=10, KTA2=10, Te1=200, Te2=200

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CONTAINS

SUBROUTINE FILE_NAMES()
	
	inFile 		= 	"in.txt"
	geomFile 	= 	"geom.txt"
	outFile 	=	"out.txt"
	ext 		= 	" .xls"
	
END SUBROUTINE FILE_NAMES
	
	
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Инициализация данных
SUBROUTINE DATA_INIT()

	INTEGER :: ti
	
	IF (OS_KEY == 1) THEN
		delimiter = '\'
	ELSE
		delimiter = '/'
	END IF

	CALL READ_DATA()

	WRITE(*,*) "Array allocation..."
	WRITE(*,11) "T(0:", Mi, ", 0:", Mj, ", 0:", Mk, ", ", prevCnt, ":", layCnt, ")"
	11 FORMAT (a,i3,a,i3,a,i3,a,i3,a,i3,a)
	WRITE(*,*) "Overall ", Mi*Mj*Mk*(layCnt-prevCnt), " nodes"

	ALLOCATE( T(0:Mi , 0:Mj , 0:Mk, prevCnt:layCnt)	)
	ALLOCATE( alpha1(1:Mi), beta1(1:Mi) )
	ALLOCATE( alpha2(1:Mj), beta2(1:Mj) )
	ALLOCATE( alpha3(1:Mk), beta3(1:Mk) )

	! задание первичного поля температур:
	DO ti = 0, layCnt
		DO i = 0, Mi
			DO j = 0, Mj
				DO k = 0, Mk
					T(i,j,k,ti) = T0
				END DO
			END DO
		END DO
	END DO
		
END SUBROUTINE DATA_INIT

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Чтение данных из файла
SUBROUTINE READ_DATA()

	INTEGER, PARAMETER 		:: 		geom = 29, input = 28
	CHARACTER*150 		 	::		command

	WRITE(*,*) "Reading data..."

	IF (OS_KEY == 1) THEN
		! удаляем директорию
		command='rmdir /s /q time'
		CALL system(command)
	ELSE
		! удаляем директорию
		command='rm -rf time'
		CALL system(command)
	END IF

	! создаём заново
	command='mkdir time'
	CALL system(command)

	! открыть общий файл входных данных
	OPEN(UNIT=input , FILE=inFile)

	READ(input,*)	Tmax
	READ(input,*)	tprint
	READ(input,*)	hx,hy,hz
	READ(input,*)	tau
	READ(input,*)	T0
	READ(input,*)	timeFold

	! вспомогательные вычисления
	hx2 = hx**2
	hy2 = hy**2
	hz2 = hz**2

	CLOSE(input)

	WRITE(*,*) "Reading geometry..."

	! открыть файл геометрии
	OPEN(UNIT=geom , FILE=geomFile)

	READ(geom,*)	lx, ly, lz	! габариты
	READ(geom,*)	delta1		! толщина верхней пластины
	
	READ(geom,*)	AT1			! значение температуропроводности для "верхней" пластины
	READ(geom,*)	LT1			! значение температуропроводности для "верхней" пластины
	READ(geom,*)	AT2			! значение температуропроводности для "нижней" пластины
	READ(geom,*)	LT2			! значение температуропроводности для "нижней" пластины
								
	READ(geom,*)	pCnt		! количество "мощностей"

	! расчёт количества точек
	Mi = lx/hx
	Mj = ly/hy
	Mk = lz/hz
	Mt = Tmax / tau + 1

	WRITE(output, *) "Количество точек: Mi = ", Mi, &
	&	"Mj = ", Mj, "Mk = ", Mj,  "Mt = ", Mt

	IF (pCnt > 0) THEN
		ALLOCATE(powers(0:pCnt-1))
		DO i = 0, pCnt-1
			! - читаем значения в формате [x0, y0, lx, ly, P]
			READ(geom,*)	powers(i)%x0, powers(i)%y0, powers(i)%lx, powers(i)%ly, powers(i)%Pow
		END DO
	END IF	

	CLOSE(geom)

END SUBROUTINE READ_DATA

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Завершение работы
SUBROUTINE DATA_FINALIZE()

	WRITE(*,*) "Finalization..."
	
	DEALLOCATE(T)
	DEALLOCATE(alpha2, beta2)
	DEALLOCATE(alpha3, beta3)
	DEALLOCATE(powers)
	
END SUBROUTINE DATA_FINALIZE




END MODULE

