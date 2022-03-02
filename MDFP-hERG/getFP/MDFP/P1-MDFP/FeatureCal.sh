echo 0 0 | gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact
echo 2 2 | gmx rms -s em.tpr -f md_0_10_center.xtc -tu ns -o rmsd.xvg
echo Bond | gmx energy -f md_0_10.edr -o enrg1.xvg
echo Angle | gmx energy -f md_0_10.edr -o enrg2.xvg
echo Proper-Dih. | gmx energy -f md_0_10.edr -o enrg3.xvg
echo LJ-14 | gmx energy -f md_0_10.edr -o enrg4.xvg
echo Coulomb-14 | gmx energy -f md_0_10.edr -o enrg5.xvg
echo LJ-\(SR\) | gmx energy -f md_0_10.edr -o enrg6.xvg
echo Coulomb-\(SR\) | gmx energy -f md_0_10.edr -o enrg7.xvg
echo Coul.-recip. | gmx energy -f md_0_10.edr -o enrg8.xvg
echo Potential | gmx energy -f md_0_10.edr -o enrg9.xvg
echo Kinetic-En. | gmx energy -f md_0_10.edr -o enrg10.xvg
echo Total-Energy | gmx energy -f md_0_10.edr -o enrg11.xvg
echo Conserved-En. | gmx energy -f md_0_10.edr -o enrg12.xvg
echo Temperature | gmx energy -f md_0_10.edr -o enrg13.xvg
echo Pressure | gmx energy -f md_0_10.edr -o enrg14.xvg
echo Constr.-rmsd | gmx energy -f md_0_10.edr -o enrg15.xvg
echo Box-X | gmx energy -f md_0_10.edr -o enrg16.xvg
echo Box-Y | gmx energy -f md_0_10.edr -o enrg17.xvg
echo Box-Z | gmx energy -f md_0_10.edr -o enrg18.xvg
echo Volume | gmx energy -f md_0_10.edr -o enrg19.xvg
echo Density | gmx energy -f md_0_10.edr -o enrg20.xvg
echo pV | gmx energy -f md_0_10.edr -o enrg21.xvg
echo Enthalpy | gmx energy -f md_0_10.edr -o enrg22.xvg
echo Vir-XX | gmx energy -f md_0_10.edr -o enrg23.xvg
echo Vir-XY | gmx energy -f md_0_10.edr -o enrg24.xvg
echo Vir-XZ | gmx energy -f md_0_10.edr -o enrg25.xvg
echo Vir-YX | gmx energy -f md_0_10.edr -o enrg26.xvg
echo Vir-YY | gmx energy -f md_0_10.edr -o enrg27.xvg
echo Vir-YZ | gmx energy -f md_0_10.edr -o enrg28.xvg
echo Vir-ZX | gmx energy -f md_0_10.edr -o enrg29.xvg
echo Vir-ZY | gmx energy -f md_0_10.edr -o enrg30.xvg
echo Vir-ZZ | gmx energy -f md_0_10.edr -o enrg31.xvg
echo Pres-XX | gmx energy -f md_0_10.edr -o enrg32.xvg
echo Pres-XY | gmx energy -f md_0_10.edr -o enrg33.xvg
echo Pres-XZ | gmx energy -f md_0_10.edr -o enrg34.xvg
echo Pres-YX | gmx energy -f md_0_10.edr -o enrg35.xvg
echo Pres-YY | gmx energy -f md_0_10.edr -o enrg36.xvg
echo Pres-YZ | gmx energy -f md_0_10.edr -o enrg37.xvg
echo Pres-ZX | gmx energy -f md_0_10.edr -o enrg38.xvg
echo Pres-ZY | gmx energy -f md_0_10.edr -o enrg39.xvg
echo Pres-ZZ | gmx energy -f md_0_10.edr -o enrg40.xvg
echo '#Surf*SurfTen' | gmx energy -f md_0_10.edr -o enrg41.xvg
echo Box-Vel-XX | gmx energy -f md_0_10.edr -o enrg42.xvg
echo Box-Vel-YY | gmx energy -f md_0_10.edr -o enrg43.xvg
echo Box-Vel-ZZ | gmx energy -f md_0_10.edr -o enrg44.xvg
echo T-System | gmx energy -f md_0_10.edr -o enrg45.xvg
echo Lamb-System | gmx energy -f md_0_10.edr -o enrg46.xvg
echo 2 3 | gmx hbond -s md_0_10.tpr -f md_0_10.xtc -num fws_hnum.xvg
gmx gyrate -s md_0_10.tpr -f md_0_10.xtc -o fws-gyrate.xvg << EOF
2
EOF
gmx sasa -s md_0_10.tpr -f md_0_10.xtc -o area.xvg -or resarea.xvg -oa atomarea.xvg << EOF
2
EOF