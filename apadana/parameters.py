def interfaceemami():
    """this function just returns emami coeffs for interface force field."""
    potential = '''


    ;;[ defaults ]
    ; nbfunc	comb-rule	gen-pairs	fudgeLJ	fudgeQQ
    ;;   1	2	yes	1.000000	1.000000




    ; Interface FF interfaceEmami par   

    [ atomtypes ]
    ; name	at.num	mass	charge	ptype	sigma	epsilon	;	sigma_14	epsilon_14 ;;
        IHOY     1     1.0080      0.400     A    9.66625109182e-02    6.276000e-02 
       IOC23     8    15.9994     -0.550     A    3.09141855195e-01    2.259360e-01 
       IOC24     8    15.9994     -0.675     A    3.09141855195e-01    5.104480e-01 
       IOC25     8    15.9994     -0.655     A    0.315378146222       0.6363864      ;  ;; 3.09141855195e-01    5.104480e-01
        ISC4    14    28.0860      1.100     A    3.69722968028e-01    3.891120e-01   ;28.085500


    [ bondtypes ]
    ; i	j	func	b0	Kb
      IOC24    IHOY     1  9.450000e-02  4.142160e+05
       ISC4   IOC23     1  1.680000e-01  2.384880e+05
       ISC4   IOC24     1  1.680000e-01  2.384880e+05
       ISC4   IOC25     1  1.680000e-01  2.384880e+05


    [ pairtypes ]
    ; i	j	func	sigma1-4	epsilon1-4


    [ angletypes ]
    ; i	j	k	func	th0	Kth	s0	Kub
      IOC23    ISC4   IOC23     5  1.0950000e+02  8.3680000e+02  0.0000000e+00  0.0000000e+00
      IOC23    ISC4   IOC24     5  1.0950000e+02  8.3680000e+02  0.0000000e+00  0.0000000e+00
      IOC23    ISC4   IOC25     5  1.0950000e+02  8.3680000e+02  0.0000000e+00  0.0000000e+00
      IOC24    ISC4   IOC24     5  1.0950000e+02  8.3680000e+02  0.0000000e+00  0.0000000e+00
      IOC24    ISC4   IOC25     5  1.0950000e+02  8.3680000e+02  0.0000000e+00  0.0000000e+00
      IOC25    ISC4   IOC25     5  1.0950000e+02  8.3680000e+02  0.0000000e+00  0.0000000e+00
       ISC4   IOC23    ISC4     5  1.4900000e+02  8.3680000e+02  0.0000000e+00  0.0000000e+00
       ISC4   IOC24    IHOY     5  1.1500000e+02  4.1840000e+02  0.0000000e+00  0.0000000e+00



    [ dihedraltypes ]
    ; i	j	k	l	func	phi0	Kphi	mult
       IHOY   IOC24    ISC4   IOC24     9  0.000000e+00  0.000000e+00      3
       IHOY   IOC24    ISC4   IOC25     9  0.000000e+00  0.000000e+00      3
      IOC23    ISC4   IOC24    IHOY     9  0.000000e+00  0.000000e+00      3
      IOC24    ISC4   IOC23    ISC4     9  0.000000e+00  0.000000e+00      3
      IOC25    ISC4   IOC23    ISC4     9  0.000000e+00  0.000000e+00      3
       ISC4   IOC23    ISC4   IOC23     9  0.000000e+00  0.000000e+00      2
    '''
    return potential


def charmm_extra():
    with open('ff_charmm_connection_point.itp', 'w+') as c:
        # Write bond parameters with unit conversion and formatted spacing
        c.write("\n\n\n\n[ bondtypes ]\n")
        c.write("; i	j	func	b0	Kb\n")
        c.write(
            "    SI    CG321       1   0.18647000    132906.43      ;260829.90=*si/c       \n")  # done  0.18647000  #;
        c.write("    SI    IOC25       1   0.16980000    252713.60             \n")  # done

        # Write angle parameters with unit conversion and formatted spacing
        c.write("\n[ angletypes ]\n")
        c.write("; i	j	k	func	th0	Kth	s0	Kub\n")
        c.write(
            "    SI      CG321    CG321       5   114.660000   348.543380   0.00000000         0.00     ;;113.600000   488.272800   0.2561000       9338.69  ;;684.019090  \n")  # done; ;;113.6 488.272800  0.2561000  9338.69
        c.write(
            "    IOC25      SI    CG321       5   111.500000   376.560000   0.00000000         0.00     ;;111.500000   376.560000   0.00000000         0.00    \n")  # C C O  ;;111.500000   376.560000   0.00000000         0.00
        c.write(
            "    OSIH       SI    CG321       5   110.100000   633.457600   0.00000000         0.00     ;;110.100000   633.457600   0.00000000         0.00    \n")  # C C O  ;;110.100000   633.457600   0.00000000         0.00
        c.write(
            "    SI      CG321    HGA2        5   110.100000   221.752000   0.21790000     18853.10     ;;110.100000   221.752000   0.21790000     18853.10    \n")  # done   ;;110.100000   221.752000   0.21790000     18853.10
        c.write(
            "    ISC4    IOC25    SI          5   150.000000   285.512000   0.00000000         0.00         \n")  # done
        c.write(
            "    IOC25      SI    OSIH        5   126.000000   267.776000   0.00000000         0.00         \n")  # done

        # Write dihedral parameters with unit conversion and formatted spacing
        c.write("\n[ dihedraltypes ]\n")
        c.write("; i	j	k	l	func	phi0	Kphi	mult\n")
        c.write("    ISC4    IOC25    SI     CG321      9     0.000000     2.384880     1       \n")  # C O C C
        c.write("    ISC4    IOC25    SI     CG321      9     0.000000     1.213360     2       \n")  # C O C C
        c.write("    ISC4    IOC25    SI     CG321      9     0.000000     1.799120     3       \n")  # C O C C

        c.write("    HSIO     OSIH    SI     CG321      9     0.000000     4.727920     1       \n")  # done
        c.write("    HSIO     OSIH    SI     CG321      9     0.000000     0.585760     2       \n")  # done
        c.write("    HSIO     OSIH    SI     CG321      9     0.000000     1.004160     3       \n")  # done

        c.write("    ISC4    IOC25    SI     OSIH       9     0.000000     0.836800     3       \n")  # q

        c.write("    IOC25    SI    CG321     CG321     9   180.000000     0.669440     1       \n")  # O C C C
        c.write("    IOC25    SI    CG321     CG321     9     0.000000     1.631760     2       \n")  # O C C C

        c.write("    IOC25    SI    CG321     HGA2      9     0.000000     0.794960     3       \n")  # O C C H

        c.write("    IOC25      SI     OSIH     HSIO    9     0.000000     1.255200     3       \n")  # q

        c.write("    IOC23    ISC4    IOC25     SI      9     0.000000     0.000000     3       \n")  # Interface

        c.write("    IOC24    ISC4    IOC25     SI      9     0.000000     0.502080     5       \n")  # q

        c.write("    IOC25    ISC4    IOC25     SI      9     0.000000     0.753120     5       \n")  # q

        # c.write("     HGA2   CG321   CG321   CG331      9     0.000000   7.949600e-01   3       \n") # done  back

        c.write("    OSIH    SI    CG321     CG321      9     0.000000     0.815880     3       \n")  # done

        c.write("    OSIH     SI     CG321    HGA2      9     0.000000     0.815880     3       \n")  # done

        c.write("    OSIH     SI     OSIH     HSIO      9     0.000000     1.464400     3       \n")  # q

        c.write("    SI    CG321     CG321    CG321     9     0.000000     0.269868     2       \n")  # done
        c.write("    SI    CG321     CG321    CG321     9   180.000000     0.626554     3       \n")  # done
        c.write("    SI    CG321     CG321    CG321     9     0.000000     0.395723     4       \n")  # done
        c.write("    SI    CG321     CG321    CG321     9     0.000000     0.470742     5       \n")  # done

        c.write("    SI    CG321    CG321     HGA2      9     0.000000     0.815800     3       \n")  # done


def charmm_cccc():
    charmmcc = '''
    ;;[ defaults ]
    ; nbfunc	comb-rule	gen-pairs	fudgeLJ	fudgeQQ
    ;;1	2	yes	1.0	1.0

    [ atomtypes ]
    ; name	at.num	mass	charge	ptype	sigma	epsilon	;	sigma_14	epsilon_14
        CG321     6    12.0110     -0.180     A    3.58141284692e-01    2.343040e-01 ;   3.38541512893e-01    4.184000e-02 
        CG331     6    12.0110     -0.270     A    3.63486677001e-01    3.263520e-01 ;   3.38541512893e-01    4.184000e-02 
        HGA2     1     1.0080      0.090     A    2.38760856462e-01    1.171520e-01 
        HGA3     1     1.0080      0.090     A    2.38760856462e-01    1.004160e-01 
    [ bondtypes ]
    ; i	j	func	b0	Kb
       CG321    CG321     1  1.530000e-01  1.861880e+05
       CG321    CG331     1  1.528000e-01  1.861880e+05
       CG321    HGA2     1  1.111000e-01  2.585712e+05
       CG331    CG331     1  1.530000e-01  1.861880e+05
       CG331    HGA3     1  1.111000e-01  2.694496e+05
    [ pairtypes ]
    ; i	j	func	sigma1-4	epsilon1-4
       CG321    CG321     1  3.38541512893e-01  4.18400000000e-02 
       CG321    CG331     1  3.38541512893e-01  4.18400000000e-02 
       CG321    HGA2     1  2.88651184677e-01  7.00117110204e-02 
       CG321    HGA3     1  2.88651184677e-01  6.48182492821e-02 
       CG331    CG331     1  3.38541512893e-01  4.18400000000e-02 
       CG331    HGA2     1  2.88651184677e-01  7.00117110204e-02 
       CG331    HGA3     1  2.88651184677e-01  6.48182492821e-02 
    [ angletypes ]
    ; i	j	k	func	th0	cth	S0	Kub
       CG321    CG321    CG321     5  1.1360000e+02  4.8827280e+02  2.5610000e-01  9.3386880e+03
       CG321    CG321    CG331     5  1.1500000e+02  4.8534400e+02  2.5610000e-01  6.6944000e+03
       HGA2    CG321    CG321     5  1.1010000e+02  2.2175200e+02  2.1790000e-01  1.8853104e+04
       HGA2    CG321    CG331     5  1.1010000e+02  2.8953280e+02  2.1790000e-01  1.8853104e+04
       HGA2    CG321    HGA2     5  1.0900000e+02  2.9706400e+02  1.8020000e-01  4.5187200e+03
       HGA3    CG331    CG321     5  1.1010000e+02  2.8953280e+02  2.1790000e-01  1.8853104e+04
       HGA3    CG331    CG331     5  1.1010000e+02  3.1380000e+02  2.1790000e-01  1.8853104e+04
       HGA3    CG331    HGA3     5  1.0840000e+02  2.9706400e+02  1.8020000e-01  4.5187200e+03
    [ dihedraltypes ]
    ; i	j	k	l	func	phi0	cp	mult
       CG321    CG321    CG321    CG321     9  0.000000e+00  3.096160e-01      4
       CG321    CG321    CG321    CG321     9  0.000000e+00  4.058480e-01      5
       CG321    CG321    CG321    CG321     9  0.000000e+00  4.225840e-01      2
       CG321    CG321    CG321    CG321     9  1.800000e+02  5.941280e-01      3
       CG321    CG321    CG321    CG331     9  0.000000e+00  4.393200e-01      4
       CG321    CG321    CG321    CG331     9  0.000000e+00  6.778080e-01      2
       CG321    CG321    CG321    CG331     9  0.000000e+00  7.405680e-01      5
       CG321    CG321    CG321    CG331     9  1.800000e+02  1.966480e-01      3
       CG331    CG321    CG321    CG331     9  0.000000e+00  1.464400e-01      5
       CG331    CG321    CG321    CG331     9  0.000000e+00  2.510400e-01      2
        HGA2    CG321    CG321     HGA2     9  0.000000e+00  7.949600e-01      3
        HGA2    CG321    CG331     HGA3     9  0.000000e+00  6.694400e-01      3
       CG321    CG321    CG321     HGA2     9  0.000000e+00  6.380600e-01      3
       CG321    CG321    CG331     HGA3     9  0.000000e+00  6.694400e-01      3
        HGA2    CG321    CG321    CG331     9  0.000000e+00  7.949600e-01      3

          '''
    return charmmcc
