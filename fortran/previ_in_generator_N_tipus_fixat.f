      PROGRAM previGeneratorNtipusFixat
C       PROGRAMA PER GENERAR EL PREVI.IN
C           PEL 2-TYPE PARTICLE

C Simplement afegim un index extra a les posicions
C Que tenim guardades del previ.in de 1-Type

       IMPLICIT REAL*8 (A-H, O-Z)

CC{{{ Introduim les variables necesaries
       INTEGER dmnp, dmnpop
       PARAMETER (dmnp = 500, dmnpop = 100)
       DIMENSION xaux (2, dmnp, dmnpop)
       DIMENSION elocal(dmnpop)
CC}}}
CC{{{ Fixem la SEED del ran1
      kkk = 4532
CC}}}
c
CC{{{ Llegim les dades del fitxer previ.in del 1-Type
       OPEN(1, file = 'previ_generator.in',status = 'old')
         READ(1,*) npop
         READ(1,*) np
         READ(1,*) N1max
         READ(1,*) N2max
         DO ipop = 1, npop
            READ(1,*) elocal (ipop)
              DO ip = 1, np
                 READ(1,*) xaux(1,ip,ipop), xaux(2,ip,ipop)
              END DO
         END DO
       CLOSE(1)
CC}}}
c
CC{{{ Afegim l'indicador de tipus de particula en el fitxer previ.in i
C     asignem els tipus a cada particula
       N1 = 0
       N2 = 0
       !N1max = 100
       !N2max = 100
       !np = 200
       OPEN(2,file='previ_initial.in',status = 'unknown')
       WRITE(2,*) npop
       WRITE(2,*) np
       WRITE(2,*) N1max
       WRITE(2,*) N2max
       DO ipop = 1, npop
          WRITE(2,*) elocal (ipop)
          !DO ip = 1, np
            DO i = 1, np
               r1 = ran1(kkk)
               r2 = ran1(kkk)
               r = ran1(kkk)
               IF ((r .LT. 0.5d0) .AND. (N1 .LT. N1max)) THEN
                   xaux(1,i,ipop) = -5.d0 + 4.d0*r1
                   xaux(2,i,ipop) = -5.d0 + 10.d0*r2
                   WRITE(2,*) xaux(1,i,ipop), xaux(2,i,ipop), 1
                   N1 = N1 + 1
               ELSEIF ((r .LT. 0.5d0) .AND. (N1 .EQ. N1max)) THEN
                   xaux(1,i,ipop) = 1.d0 + 4.d0*r1
                   xaux(2,i,ipop) = -5.d0 + 10.d0*r2
                   WRITE(2,*) xaux(1,i,ipop), xaux(2,i,ipop), 2
                   N2 = N2 + 1
               ELSEIF ((r .GT. 0.5d0) .AND. (N2 .LT. N2max)) THEN
                   xaux(1,i,ipop) = 1.d0 + 4.d0*r1
                   xaux(2,i,ipop) = -5.d0 + 10.d0*r2
                   WRITE(2,*) xaux(1,i,ipop), xaux(2,i,ipop), 2
                   N2 = N2 + 1
               ELSE ! ((r .GT. 0.5d0) .AND. (N2 .EQ. N2max))
                   xaux(1,i,ipop) = -5.d0 + 4.d0*r1
                   xaux(2,i,ipop) = -5.d0 + 10.d0*r2
                   WRITE(2,*) xaux(1,i,ipop), xaux(2,i,ipop), 1
                   N1 = N1 + 1
               END IF
            END DO
          END DO
        !END DO
        CLOSE(2)

CC}}}


      END PROGRAM previGeneratorNtipusFixat

c
c
c
CC{{{ ran1: Per generar random deviates U[0,1]
        function ran1(idum)
c
        implicit real*8 (a-h,o-z)
c
        SAVE
c
        dimension r(97)
        parameter ( m1=259200, ia1=7141, ic1=54773, rm1=1./m1)
        parameter ( m2=134456, ia2=8121, ic2=28411, rm2=1./m2)
        parameter ( m3=243000, ia3=4561, ic3=51349)
c
        data iff / 0 /
c
        if(idum.lt.0 .or. iff.eq.0) then
          iff = 1
          ix1 = mod ( ic1 - idum , m1 )
          ix1 = mod ( ia1 * ix1 + ic1 , m1 )
          ix2 = mod ( ix1, m2 )
          ix1 = mod ( ia1 * ix1 + ic1 , m1 )
          ix3 = mod ( ix1, m3 )
          do j = 1, 97
             ix1 = mod ( ia1 * ix1 + ic1 , m1 )
             ix2 = mod ( ia2 * ix2 + ic2 , m2 )
             r(j) = ( float(ix1) + float(ix2) * rm2 ) * rm1
          end do
        end if
c
c       ccccccccccccccccccccccc   except when initializing, this is
c                                 were we start.
        ix1 = mod ( ia1 * ix1 + ic1 , m1 )
        ix2 = mod ( ia2 * ix2 + ic2 , m2 )
        ix3 = mod ( ia3 * ix3 + ic3 , m3 )
        j = 1 + ( 97 * ix3 ) / m3
c
c
        if ( j.gt.97 .or. j.lt.1 ) write(6,*) ' AAAUUGGHH !!!'
c
        ran1 = r(j)
        r(j) = ( float(ix1) + float(ix2) * rm2 ) * rm1
c------------------------------------------------------
c       do j = 1 , 97
c       write(6,*) ' j, r(j) = ', j, r(j)
c       end do
c------------------------------------------------------
c       write(6,*)
c
        return
        end
CC}}}
