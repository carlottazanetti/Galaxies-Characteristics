module dati 
    implicit none 
    real*8,parameter:: H0=70.d0*3.24e-20 , gamma0l=0.7, gamma0m=0.3,pi=acos(-1.d0)
    integer,parameter:: npassi=1000
    save 
end module dati 
program esame
    use dati 
    implicit none
    real*8,allocatable::lambda_sed(:), lambda_gal(:),f_gal(:),err_f(:),f_sed(:), sed_finale(:), lambdased_finale(:), &
    & transmission(:),lambda_filter(:),derivate_filter(:),matr1(:,:),matr2(:,:),sed_filter(:), &
    & integrandaLT(:), derivate_integrandaLT(:), derivate_sed(:), derivate(:), f_sim(:)
    integer:: nd_sed , i , j , nd_gal , k , ndsed_finale, nd_filter, jran 
    real*8::  somma, best_chi, best_norm, z , chi_min(280), chi_grafico(280), z_grafico(280), &
    & result , redshift, chi_finale, norm_finale,z_finale,D, integrale1 ,L, &
    & integrale3, tLB, t_z, integrale4, a, b , integrale2,integraleT, integraleLT, flux_limit, &
    & flux_filter(500), D_filter(500), z_filter(500), lambda_eff, integrale5,L_sim(100),z_sim(100), &
    & rms_L, media_L, rms_z, media_z
    character (len=15):: names_sed(5),sed_nome 
    character (len=13) names_gal(4)
    character (len=14) names_saved(20)
    real*8,external:: E , E2 , E3

    !Leggo immediatamente i dati del filtro per usarli nel ciclo successivo (così non devo ripetere il calcolo per ogni ciclo)
    open(13,file='pacs_transmission_blue.txt')
    nd_filter=-1 !così i dati letti sono nd+1
    read(13,*)
    !leggo gli N+1 dati del filtro
    do 
        read(13,*,end=1002)
        nd_filter=nd_filter+1
    end do 
    1002 continue 
    rewind(13)
    allocate(transmission(nd_filter+1),lambda_filter(nd_filter+1),derivate_filter(nd_filter+1),integrandaLT(nd_filter+1))
    allocate(matr1(nd_filter+1,nd_filter+1),matr2(nd_filter+1,nd_filter+1),sed_filter(nd_filter+1))
    allocate(derivate_integrandaLT(nd_filter+1))
    do j=1,nd_filter +1
        read(13,*) lambda_filter(j), transmission(j) !lambda in micron 
    end do
    close(13)
    transmission=log10(transmission) !metto in logaritmo
    lambda_filter=log10(lambda_filter)
    !Faccio l'integrale della trasmissione dato che non necessita di essere ripetuto per ogni ciclo
    call sorting(lambda_filter,transmission,nd_filter+1)!ordino i valori di lambda in ordine crescente
    call matrice(lambda_filter,transmission,nd_filter,matr1,matr2) !matrice per fare la spline cubica
    call gauss(matr1,matr2,nd_filter+1,derivate_filter) !trovo le derivate seconde 
    a=14.d0 !Estremi delle lunghezze d'onda coinvolte dal filtro in micron
    b=140.d0
    call integrale(a,b,nd_filter,transmission,lambda_filter,derivate_filter,0.d0,integraleT)
    !Calcolo anche lambda effettivo per le conversioni
    call integralelambda(a,b,nd_filter,transmission,lambda_filter,derivate_filter,integrale5)!subroutine che interpola T nei punti da 14 a 140 micron e calcola l'integrale di lambda*T(lambda)
    lambda_eff=integrale5/integraleT
 
    names_sed(1)='S0_sed_norm.txt' !vettori con i nomi che userò nel ciclo
    names_sed(2)='Sb_sed_norm.txt'
    names_sed(3)='Sc_sed_norm.txt'
    names_sed(4)='Sd_sed_norm.txt'
    names_sed(5)='El_sed_norm.txt'
    names_gal(1)='galaxy_01.txt'
    names_gal(2)='galaxy_02.txt'
    names_gal(3)='galaxy_03.txt'
    names_gal(4)='galaxy_04.txt'
    names_saved(1)='Galassia_1.txt'
    names_saved(5)='Sed_Gal1.txt'
    names_saved(2)='Galassia_2.txt'
    names_saved(6)='Sed_Gal2.txt'
    names_saved(3)='Galassia_3.txt'
    names_saved(7)='Sed_Gal3.txt'
    names_saved(4)='Galassia_4.txt'
    names_saved(8)='Sed_Gal4.txt'
    names_saved(9)='Chi_Gal1.txt'
    names_saved(10)='Chi_Gal2.txt'
    names_saved(11)='Chi_Gal3.txt'
    names_saved(12)='Chi_Gal4.txt'
    names_saved(13)='Filtro_1.txt'
    names_saved(14)='Filtro_2.txt'
    names_saved(15)='Filtro_3.txt'
    names_saved(16)='Filtro_4.txt'
    names_saved(17)='ZL_sim_1.txt'
    names_saved(18)='ZL_sim_2.txt'
    names_saved(19)='ZL_sim_3.txt'
    names_saved(20)='ZL_sim_4.txt'

    do i=1,4
        open(11,file=names_gal(i))
        nd_gal=-1
        read(11,*)
        !leggo gli N+1 dati della galassia
        do 
            read(11,*,end=1000)
            nd_gal=nd_gal+1
        end do 
        1000 continue 
        rewind(11)
        if (allocated(lambda_gal)) deallocate(lambda_gal)
        if (allocated(f_gal)) deallocate(f_gal)
        if (allocated(err_f)) deallocate(err_f)
        if (allocated(f_sim)) deallocate(f_sim) !vettore del flusso simulato per il facoltativo
        allocate(lambda_gal(nd_gal+1),f_gal(nd_gal+1),err_f(nd_gal+1),f_sim(nd_gal+1))
        read(11,*)
        do j=1,nd_gal +1
            read(11,*) lambda_gal(j) , f_gal(j), err_f(j) !in micron e flusso per unità di frequenza in microjanksy
        end do
        close(11)
        print*, 'Sto analizzando la galassia ', names_gal(i)
        lambda_gal=lambda_gal*10.**4 !da micron a angstrom 
        f_gal=f_gal *10.d0**(-29)*(((2.99e10)*10.**8)/((lambda_gal)**2))!da Ffrequenza in Flambda (erg cm^-2 s^-1 angstrom^-1)
        err_f=err_f*10.d0**(-29)*((2.99e10)*10.**8)/((lambda_gal)**2) !conveto anche l'errore sul flusso
        !trasformo in log (non trasformo lambda perché lo trasformerò in seguito nella subroutine del chi quadro dopo averlo convertito in lambda emesso)
        err_f=(0.43*err_f)/f_gal !propagazione dell'errore in log10
        f_gal=log10(f_gal) !metto il flusso in logaritmo
        chi_finale=10.d0**8 !valore arbitrario per entrare del ciclo del chi quadro
        print*,'Sto analizzando le sed per questa galassia'
        do j=1,5
            open(12,file=names_sed(j))
            nd_sed=-1
            !leggo gli n+1 dati della sed
            do 
                read(12,*,end=1001)
                nd_sed=nd_sed+1
            end do
            1001 continue 
            rewind(12)
            if (allocated(lambda_sed)) deallocate(lambda_sed)
            if (allocated(f_sed)) deallocate(f_sed)           
            if(allocated(derivate)) deallocate(derivate)
            allocate(f_sed(nd_sed+1),lambda_sed(nd_sed+1),derivate(nd_sed+1))
            do k=1,nd_sed+1
                read(12,*) lambda_sed(k),f_sed(k) !in angstrom e flusso per unità di lambda
            end do 
            close(12)
            lambda_sed=log10(lambda_sed)!metto in log10
            f_sed=log10(f_sed)
            call chiquadro(nd_gal,lambda_gal,f_gal,nd_sed,lambda_sed,f_sed,err_f,best_chi,redshift,best_norm,derivate,chi_min) !ottengo il chi quadro, il redshift e la normalizzazione migliore per ogni sed 
            if (j==1) print*, 'Sto facendo il chi quadro'  !printo la frase solo una volta    
            if (best_chi<chi_finale) then 
                if (allocated(derivate_sed)) deallocate(derivate_sed)
                if (allocated(sed_finale)) deallocate(sed_finale)
                if (allocated(lambdased_finale)) deallocate(lambdased_finale)
                allocate(derivate_sed(nd_sed+1))
                allocate(sed_finale(nd_sed+1))
                allocate(lambdased_finale(nd_sed+1))
                derivate_sed=derivate !Le derivate mi serviranno per fare l'interpolazione in vari punti del progetto
                sed_finale=f_sed !salvo la sed migliore
                lambdased_finale=lambda_sed !salvo le lugnhezze d'onda
                ndsed_finale=nd_sed !salvo il numero di dati
                chi_finale=best_chi !salvo il chi quadro migliore
                norm_finale=best_norm  !salvo la normalizzazione migliore 
                z_finale=redshift !salvo il redshift corretto
                sed_nome=names_sed(j) !salvo il nome della sed giusta
                chi_grafico=chi_min !salvo il vettore del chi quadro della sed migliore per il grafico
            end if
        end do
        print*, 'la galassia ', names_gal(i),' corrisponde a ', sed_nome,' con un redshift di',&
        & z_finale,'e una normalizzazione di',norm_finale 
        call save_results2(names_saved(i),log10(lambda_gal/(1.d0+z_finale)),f_gal,err_f,nd_gal+1) !Salvo i dati della galassia 
        call save_results(names_saved(i+4),lambdased_finale,sed_finale-norm_finale,ndsed_finale+1) !salvo i dati della sed migliore
        do j=1,280 !mi salvo un vettore con tutti i valori del redshift per il grafico
            z=j*0.01
            z_grafico(j)=z
        end do
        call save_results(names_saved(i+8),z_grafico,chi_grafico,280)!grafico chi quadro in funzine di z

        !Distanza luminosa 
        call trapezoide(E,0.d0,z_finale,1000,integrale1) !integrale di E^-1 per la distanza
        D=(9.69e-15/H0)*(1.d0+z_finale)*integrale1 !in Mpc
        print*,'la sua distanza luminosa è ', D, ' Mpc'

        !Calcolo luminosità
        a=8.d0*10.d0**4 !Estremi dell'integrale da micron a angstrom 
        b=1000.d0*10.d0**4
        call integrale(a,b,ndsed_finale,sed_finale,lambdased_finale,derivate_sed,norm_finale,integrale2) !integrale di F in dlambda
        D=D*3.086e24 !da Megaparsec a cm
        L=4.d0*pi*D**2*(1.d0+z_finale)*integrale2 !in erg/cm
        print*, 'la sua luminosità è', L,'erg/s.'

        !ETA' DELL'UNIVERSO AL MOMENTO DELL'EMISSIONE DELLA GALASSIA E LOOK BACK TIME
        !Calcolo t(z)
        call trapezoide(E3,1.d-14,1/z_finale,1000,integrale3) !integrale che va da 1/infinito a 1/z di E^-1*(1+z)^-1 con il cambio di variabile 1/t
        t_z=(1.d0/H0)*integrale3!in secondi
        t_z=t_z*3.171e-8 !in anni
        !Calcolo il look back time 
        call trapezoide(E2,0.d0,z_finale,1000,integrale4) !integrale di E^-1*(1+z)^-1 da 0 a z
        tLB=(1.d0/H0)*integrale4 !in secondi
        tLB=tLB*3.171e-8 !in anni
        print*, "L'età dell'universo al momento dell'emissione di questa galassia era", t_z, &
        &'anni, e il look back time risulta', tLB, 'anni.'

        !SECONDA PARTE 
        lambdased_finale=10.d0**lambdased_finale
        lambdased_finale=0.0001*lambdased_finale !da angstrom in micron
        lambdased_finale=log10(lambdased_finale)
        sed_finale=10.d0**sed_finale
        sed_finale=sed_finale*10.d0**4!Da erg cm^-2 s^-1 angstrom^-1 a erg cm^-2 s^-1 micron^-1
        sed_finale=log10(sed_finale)
        do j=1,nd_filter+1
            call interpolazione(lambda_filter(j),derivate_sed,sed_finale,ndsed_finale,result,lambdased_finale) !interpolo la sed nelle lunghezze d'onda del filtro
            result=result-norm_finale !tolgo la normalizzazione
            sed_filter(j)=result
        end do
        sed_filter=10.d0**sed_filter
        transmission=10.d0**transmission
        integrandaLT=sed_filter*4.d0*pi*D**2*(1.d0+z_finale)*transmission
        sed_filter=log10(sed_filter)!rimetto tutto in log
        transmission=log10(transmission)
        integrandaLT=log10(integrandaLT)
        a=14.d0 !Estremi dell'integrale in micron
        b=140.d0
        call sorting(lambda_filter,integrandaLT,nd_filter+1)!faccio la spline 
        call matrice(lambda_filter,integrandaLT,nd_filter,matr1,matr2)
        call gauss(matr1,matr2,nd_filter+1,derivate_integrandaLT) 
        call integrale(a,b,nd_filter,integrandaLT,lambda_filter,derivate_integrandaLT,0.d0,integraleLT)
        somma=0.d0
        do j=1,500 !valori del redshift da 0.01 a 5
            z=0.01*j
            call trapezoide(E,0.d0,z,1000,integrale1) 
            D_filter(j)=(3.e10/H0)*(1.d0+z)*integrale1 !in cm
            flux_filter(j)=integraleLT/(4.d0*pi*D_filter(j)**2*(1.d0+z)*integraleT) !flusso lambda attraverso filtro per ogni z
            flux_filter(j)=(flux_filter(j)*lambda_eff**2)/(3.e8*1000000.d0) !conversione da flusso per unità de lambda a flusso per unità di frequenza
            D_filter(j)=D_filter(j)*3.086e-25 !Da cm a Mpc
            z_filter(j)=z !Salvo i dati in vettori per i grafici
            flux_limit=1.e-26 !Flusso liite da mJy a erg/s/cm^2/Hz
            if (flux_filter(j)<=flux_limit .and. somma==0.d0) then 
                print*, 'La distanza massima dalla quale questo oggetto è osservabile è', D_filter(j-1),&
                & 'Mpc, che corrisponde ad un redshift di', z_filter(j-1),&
                &'e ad un flusso per unità di frequenza di',flux_filter(j-1) ,'erg/s/cm^2/Hz'
                somma=1.d0 !stampo il risultato solo un volta quando supero il flusso limite
            end if
        end do
        call save_results(names_saved(i+12),z_filter,log10(flux_filter),500) !salvo i dati per il grafico del flusso attraverso il filtro

        !FACOLTATIVO 
        print*, 'Facoltativo: set di 100 simulazioni del flusso delle galassie'
        sed_finale=10.d0**sed_finale
        sed_finale=sed_finale*10.d0**(-4) !rimetto la sed in erg cm^-2 s^-1 angstrom^-1
        sed_finale=log10(sed_finale)
        lambdased_finale=10.d0**lambdased_finale
        lambdased_finale=lambdased_finale*10.d0**4 !rimetto in angstrom
        lambdased_finale=log10(lambdased_finale)
        jran=4858 !valore del seme iniziale
        do j=1,100 !numero di simulazioni
            f_gal=10.d0**f_gal !tolgo il logaritmo
            err_f=err_f*f_gal/0.43
            f_sim=0. 
            call montecarlo (jran,nd_gal+1,f_gal,err_f,f_sim) !calcolo i vettore dei flussi simulati
            err_f=(0.43*err_f)/f_gal
            f_gal=log10(f_gal)
            f_sim=log10(f_sim) !rimetto il log
            if (j==25) print*, 'Sto cercando il redshift per ognuna delle 100 simulazioni '
            call chiquadrosim(nd_gal,lambda_gal,f_sim,ndsed_finale,lambdased_finale,sed_finale,err_f,&
            & derivate_sed,redshift,norm_finale) 
            z_sim(j)=redshift !vettore con tutte le simulazioni del redshift
            call trapezoide(E,0.d0,redshift,1000,integrale1) !integrale di E^-1
            D=(9.69e-15/H0)*(1.d0+redshift)*integrale1 !in Mpc
            a=8.d0*10.d0**4 !Estremi dell'integrale da micron a angstrom 
            b=1000.d0*10.d0**4
            if (j==75) print*, 'Sto calcolando la luminosità per ognuna delle 100 simulazioni'
            call integrale(a,b,ndsed_finale,sed_finale,lambdased_finale,derivate_sed,norm_finale,integrale2) !integrale di F in dlambda
            D=D*3.086e24 !da Megaparsec a cm
            L=4.d0*pi*D**2*(1.d0+redshift)*integrale2
            L_sim(j)=L !vettore con le luminosità simulate 
        end do
        call media_varianza(z_sim, 100,media_z, rms_z) 
        print*, "La media dei valori del redshift simulati è ",media_z, 'con un errore gaussiano di', rms_z 
        call media_varianza(L_sim, 100,media_L, rms_L)
        print*, "La media dei valori della luminosità simulati è ", &
        & media_L, 'con un errore gaussiano di', rms_L
        call save_results(names_saved(i+16),z_sim,L_sim,100)
        print*, ' '
    end do 
    print*, 'the end'   
end program esame 

real*8 function E(z)
use dati
implicit none 
real*8:: z 
E=1.d0/(sqrt(gamma0m*(1.d0+z)**3+(1.d0-gamma0m-gamma0l)*(1.d0+z)**2+gamma0l))
end function 

real*8 function E2(z)
use dati
implicit none 
real*8:: z 
E2=1.d0/(sqrt(gamma0m*(1.d0+z)**3+(1.d0-gamma0m-gamma0l)*(1.d0+z)**2+gamma0l)*(1.d0+z))
end function 

real*8 function E3(t) 
implicit none 
real*8:: t 
real*8, external::E2 
E3=E2(1.d0/t)/(t**2)
!Ho fatto il cambio i variabile, e il dx diventa -1/t**2 dt , ma il meno si elimina con il cambio degli estremi
end function E3

subroutine sorting(col1,col2,ndati)
    implicit none
    integer, intent(in):: ndati
    real*8:: col1(ndati), col2(ndati) 
    real*8:: aux1(ndati), aux2(ndati) , min ,provv
    integer:: i, posmin , j

    aux1=col1
    aux2=col2
    do i=1,ndati-1
        min=aux1(i)
        posmin=i 
        do j=i+1,ndati 
            if(aux1(j)<min) then
                min=aux1(j)
                posmin=j
            end if
        end do
        aux1(posmin)=aux1(i)
        aux1(i)=min 
        provv=aux2(posmin)
        aux2(posmin)=aux2(i)
        aux2(i)=provv
    end do
    !Ora rimetto aux 1 e 2 nelle colonne col1 e col2
    col1=aux1
    col2=aux2

end subroutine sorting

subroutine matrice(x,f,nd,a,c) !subroutine che calcola la matrice (n+1,n+1) del sistema per trovare le derivate seconde
    implicit none 
    real*8, intent(in):: x(nd+1), f(nd+1)
    integer, intent(in):: nd
    real*8, intent(out):: a(nd+1, nd+1), c(nd+1) 
    integer::i, j 

    a=0.d0
    c=0.d0 
    do i=2,nd
        do j=1, nd+1
            if (i==j) then !diagonale principale
                a(i,j)=2*(x(i+1)-x(i-1)) 
            else if (i+1==j) then !diagonale sopra la diagonale principale
                a(i,j)=(x(i+1)-x(i))
            else if (i==j+1) then !diagonale sotto la diagonale principale
                a(i,j)=(x(i)-x(i-1))
            end if
        end do
        c(i)=6./(x(i+1)-x(i))*(f(i+1)-f(i))+(6./(x(i)-x(i-1)))*(f(i-1)-f(i)) !vettore del termine noto
    end do
    a(1,1)=1 !condizione: derivate seconde nei nodi finali sono assunte nulle
    a(nd+1,nd+1)=1

end subroutine matrice

subroutine gauss(a,c,n,derivate) !matrice che risolve il sistema algebrico e trova le derivate seconde 
    implicit none
    integer:: n, i, j, k , l
    real*8:: a(n,n), c(n), derivate(n), fakt, somma

    !Inizio la forward elimination
    do i=1,n-1 !Scelgo la variabile da eliminare
        do j=i+1,n !Scelgo la riga dove eliminare
            fakt= a(j,i)/a(i,i) !Fattore di normalizzazione
            do k=1,n !Agisco su tutte le colonne
                a(j,k)=a(j,k)-a(i,k)*fakt
            end do
            c(j)=c(j)-c(i)*fakt !Agisco sul termine noto
        end do
    end do

    !Inizio la backward sostitution
    derivate(n)=c(n)/a(n,n)
    do i=n-1,1, -1
        somma=c(i)
        do j=i+1,n 
            somma=somma-a(i,j)*derivate(j)
        end do
        derivate(i)=somma/a(i,i)
    end do

end subroutine gauss

subroutine interpolazione(x_int,derivate,f,nd,fun,estremi)
    implicit none
    integer, intent(in):: nd
    real*8:: x_int, estremi(nd+1),fun,derivate(nd+1), f(nd+1)
    integer::  i 

    if (x_int<minval(estremi) .or. x_int>maxval(estremi)) fun=0.d0  !estrapolazione
    do i=1,nd
        if (x_int>=estremi(i) .and. x_int<=estremi(i+1)) then 
            fun=(derivate(i)/(6.*(estremi(i+1)-estremi(i))))*(estremi(i+1)-x_int)**3 + &
                (derivate(i+1)/(6.*(estremi(i+1)-estremi(i))))*(x_int-estremi(i))**3 + &
                ((f(i)/(estremi(i+1)-estremi(i)))-((derivate(i)*(estremi(i+1)-estremi(i)))/6.))*(estremi(i+1)-x_int)+&
                ((f(i+1)/(estremi(i+1)-estremi(i)))-((derivate(i+1)*(estremi(i+1)-estremi(i)))/6.))*(x_int-estremi(i))  
        end if 
    end do
    
end subroutine interpolazione

subroutine chiquadro(nd_gal,lambda_gal,f_gal,nd_sed,lambda_sed,f_sed,err_f,best_chi,redshift,best_norm,derivate,chi_min)
    use dati
    implicit none 
    integer:: i, j, k 
    integer, intent(in)::  nd_gal, nd_sed
    real*8, intent(in):: lambda_gal(nd_gal+1), f_gal(nd_gal+1), err_f(nd_gal+1)
    real*8, intent(out):: best_chi, redshift, best_norm, derivate(nd_sed+1), chi_min(280)
    real*8:: lambda_em(nd_gal+1), f_sed(nd_sed+1),inter_sed(nd_gal+1), f_norm(nd_gal+1), &
    & result , diff_max, diff_min,passo, norm, z, somma , chi(npassi), diff(nd_gal+1), &
    & normalizzazione(npassi),norm_min(280),a(nd_sed+1,nd_sed+1), c(nd_sed+1), lambda_sed(nd_sed+1)

    call sorting(lambda_sed,f_sed,nd_sed+1)
    call matrice(lambda_sed,f_sed,nd_sed,a,c)
    call gauss(a,c,nd_sed+1,derivate)
    best_chi=100000.
    do i=1,280 
        z=0.01*i 
        do j=1,nd_gal+1
            lambda_em(j)=lambda_gal(j)/(1.d0+z)
            lambda_em(j)=log10(lambda_em(j))         
            call interpolazione(lambda_em(j),derivate,f_sed,nd_sed,result,lambda_sed)  !Interpolo la sed nei valori di lambda emesso 
            inter_sed(j)=result !creo un vettore con tutti i valori di f interpolati 
            diff(j)= inter_sed(j)-f_gal(j) 
        end do
        diff_max=maxval(diff)
        diff_min=minval(diff)
        passo=(diff_max-diff_min)/npassi 
        do j=1,npassi  !ciclo per ogni normalizzazione 
            norm= diff_min+j*passo 
            f_norm=f_gal+norm 
            somma=0.d0
            do k=1, nd_gal+1 !ciclo per calcolare il chi quadro per ogni normalizazione 
                somma=somma+(((f_norm(k)-inter_sed(k))**2)/((err_f(k)))**2)
            end do
            chi(j)=somma 
            normalizzazione(j)=norm 
        end do 
        call sorting(chi,normalizzazione,npassi) !ordino il chi quadro in ordine crescente
        chi_min(i)=chi(1) !chi quadro minimo per ogni z 
        norm_min(i)=normalizzazione(1) ! normalizzazione che corrisponde al chi migliore per ogni z 
        if (chi_min(i)<=best_chi) then 
            best_chi=chi_min(i) !chi migliore tra tutti i redshift
            best_norm=norm_min(i)
            redshift=z
        end if
    end do

end subroutine chiquadro

subroutine trapezoide(fun,a,b,n,result)
    implicit none
    real*8, external :: fun
    real*8:: a,b,result, delta, x1, x2
    integer::n , i 
    delta=(b-a)/n
    result=0.d0
    do i=1,n
        x1=a+(i-1)*delta
        x2=a+i*delta
        result=result+0.5*(fun(x1)+fun(x2))
    end do
    result=result*delta
end subroutine trapezoide

subroutine integrale(a,b,nd,f,lambda,derivate,normalizzazione,integrale2)
    use dati
    implicit none
    real*8:: a, b, f(nd+1), lambda(nd+1), derivate(nd+1), h , somma, f_intera, f_interb, &
    & normalizzazione , result, integrale2
    integer:: nd, j  

    h=(b-a)/10000.d0
    somma=0.d0
    call interpolazione(log10(a),derivate,f,nd,f_intera,lambda) !interpolo i due estremi
    call interpolazione(log10(b),derivate,f,nd,f_interb,lambda)   
    do j=1,10000-1 !interpolo i punti interni 
        call interpolazione(log10(a+j*h),derivate,f,nd,result,lambda) 
        if (result==0.d0) then !se cado nell'estrapolazione il valore va scartato (e non tolgo il log perché 10**0=1)
            result=0.d0
        else
            result=result-normalizzazione !tolgo la normalizzazione
            result=10.d0**result !tolgo il log 
            somma=somma+result
        end if
    end do
    if (f_intera==0.d0) then
        f_intera=0.d0
    else
        f_intera=f_intera-normalizzazione !Tolgo la normalizzazione
        f_intera=10.d0**f_intera  !tolgo il log
    end if
    if (f_interb==0.d0) then
        f_interb=0.d0
    else
        f_interb=f_interb-normalizzazione
        f_interb=10.d0**f_interb
    end if
    integrale2=(b-a)*((f_intera + 2.d0*somma+f_interb)/(2.d0*10000.d0))

end subroutine integrale

subroutine integralelambda(a,b,nd,f,lambda,derivate,integrale)
    use dati
    implicit none
    real*8:: a, b, f(nd+1), lambda(nd+1), derivate(nd+1), h , somma, f_intera, f_interb, &
    & result, integrale
    integer:: nd, j  

    call interpolazione(log10(a),derivate,f,nd,f_intera,lambda) !interpolo i due estremi
    if (f_intera==0.d0) then
        f_intera=0.d0
    else
        f_intera=10.d0**f_intera
        f_intera=f_intera*a
    end if
    call interpolazione(log10(b),derivate,f,nd,f_interb,lambda)
    if (f_interb==0.d0) then
        f_interb=0.d0
    else   
        f_interb=10.d0**f_interb
        f_interb=f_interb*b
    end if
    h=(b-a)/10000.d0
    somma=0.d0
    do j=1,10000-1 !interpolo i punti interni 
        call interpolazione(log10(a+j*h),derivate,f,nd,result,lambda)
        if (result==0.d0) then
            result=0.d0
        else
            result=10.d0**result !tolgo il log 
            result=result*(a+j*h)
            somma=somma+result
        end if
    end do
    integrale=(b-a)*((f_intera + 2.d0*somma+f_interb)/(2.d0*10000.d0))

end subroutine integralelambda

subroutine save_results(filename,x,y,nd)!subroutine che salva i risultati in 2 colonne
    integer, intent(in):: nd
    character (len=*), intent(in):: filename 
    real*8:: x(nd), y(nd)
    integer:: i 

    open(40,file=filename)
    do i=1,nd 
        write(40,'(2(1pe20.10))') x(i), y(i)
    end do
    close(40)
end subroutine save_results

subroutine save_results2(filename,x,y,err,nd)!subroutine che salva i risultati in 3 colonne
    integer, intent(in):: nd
    character (len=*), intent(in):: filename 
    real*8:: x(nd), y(nd), err(nd)
    integer:: i 

    open(40,file=filename)
    do i=1,nd 
        write(40,'(3(1pe20.10))') x(i), y(i), err(i)
    end do
    close(40)
end subroutine save_results2

subroutine montecarlo (jran,nd,f,err,f_sim)
    implicit none
    integer, parameter:: m=259200 ,c=54773 , a=7141
    integer:: nd, j , jran
    real*8:: x1,x2,y1,f(nd),err(nd),f_sim(nd)
    do j=1,nd 
        x1=0.d0
        do while (x1<=0.d0) !esco dal ciclo solo quando il valore di x1 sarà maggiore di 0 (x1 va da 0 a 1 non compreso)
            jran=mod(jran*a+c,m)
            x1=float(jran)/float(m)
        end do
        jran=mod(jran*a+c,m)
        x2=float(jran)/float(m) !altro numero casuale
        y1=sqrt(-2.d0*log(x1))*sin(2.*acos(-1.d0)*x2)
        f_sim(j)=f(j)+err(j)*y1 !valore simulato secondo una distribuzione gaussiana con media=f(j) e varianza=err(j)
    end do
end subroutine montecarlo

subroutine media_varianza(vettore, n, media, rms)
    implicit none
    integer, intent(in):: n 
    real*8, intent(in):: vettore(n) 
    real*8:: media, varianza, rms
    integer:: i 
    real*8:: somma

    somma=0.d0
    do i=1,n 
        somma=somma+vettore(i)
    end do
    media=somma/n 

    varianza=0.d0
    do i=1,n 
        varianza=varianza+(vettore(i)-media)**2
    end do
    varianza=varianza/(n-1)
    rms=sqrt(varianza)

end subroutine media_varianza

subroutine chiquadrosim(nd_gal,lambda_gal,f_gal,nd_sed,lambda_sed,f_sed,err_f,derivate,redshift,best_norm)
    use dati
    implicit none 
    integer:: i, j, k 
    integer, intent(in)::  nd_gal, nd_sed
    real*8, intent(in):: lambda_gal(nd_gal+1), f_gal(nd_gal+1), err_f(nd_gal+1), derivate(nd_sed+1)
    real*8, intent(out):: redshift, best_norm
    real*8:: lambda_em(nd_gal+1), f_sed(nd_sed+1),inter_sed(nd_gal+1), f_norm(nd_gal+1), &
    & result , diff_max, diff_min,passo, norm, z, somma , chi(npassi), diff(nd_gal+1), &
    & normalizzazione(npassi),norm_min(280),a(nd_sed+1,nd_sed+1), c(nd_sed+1), lambda_sed(nd_sed+1),&
    &  best_chi, chi_min(280)

    !Questa subroutine è uguale alla subroutine chiquadro, però non ricalcola le derivate seconde per la spline perché le ho già
    best_chi=100000.
    do i=1,280 
        z=0.01*i 
        do j=1,nd_gal+1
            lambda_em(j)=lambda_gal(j)/(1.d0+z)
            lambda_em(j)=log10(lambda_em(j))         
            call interpolazione(lambda_em(j),derivate,f_sed,nd_sed,result,lambda_sed)  
            inter_sed(j)=result  
            diff(j)= inter_sed(j)-f_gal(j) 
        end do
        diff_max=maxval(diff)
        diff_min=minval(diff)
        passo=(diff_max-diff_min)/npassi 
        do j=1,npassi   
            norm= diff_min+j*passo 
            f_norm=f_gal+norm 
            somma=0.d0
            do k=1, nd_gal+1  
                somma=somma+(((f_norm(k)-inter_sed(k))**2)/((err_f(k)))**2)
            end do
            chi(j)=somma 
            normalizzazione(j)=norm 
        end do 
        call sorting(chi,normalizzazione,npassi) 
        chi_min(i)=chi(1)
        norm_min(i)=normalizzazione(1) 
        if (chi_min(i)<=best_chi) then 
            best_chi=chi_min(i) 
            best_norm=norm_min(i)
            redshift=z
        end if
    end do

end subroutine chiquadrosim




