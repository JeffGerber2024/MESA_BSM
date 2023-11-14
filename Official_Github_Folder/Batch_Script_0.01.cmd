#! /bin/sh
#The following lines set cluster environment variables such as how many tasks to run at once (ntasks), how many cpus to use per task(cpus-per-task), how many nodes to use overall (nodes), and where to save the output and error files.
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --output=test-%j.out
#SBATCH --error=test-%j.err
#SBATCH --cpus-per-task=28

#Creation of two arrays telling which stellar masses to run as high mass and which to run as low mass
masseslow=( 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 7.0 8.0 9.0 10.0 )
masseshigh=( 12.0 14.0 16.0 18.0 20.0 )

#If the results folder does not exist, create it
if [[ ! -d 'results' ]]; then mkdir 'results'; fi

#Loop iterates over the masseslow array
for w in "${masseslow[@]}"
do 
    sleep 5
    echo "$w" #prints the current mass that is being simulated
    sed -i '$ d' inlist_all #Removes the last 4 lines of the inlist_all file (those corresponding to initial mass, helium abundance, and metallicity of the star)
    sed -i '$ d' inlist_all
    sed -i '$ d' inlist_all
    sed -i '$ d' inlist_all
    echo "initial_z = 0.01" >> inlist_all #Replaces those 4 lines with a star of mass w, metal content 0.01, and helium abundance that is solar scaled to 0.01 (Calculated by linearly scaling to the Helium abundance and metallicity that best reproduces solar observables at 4.61 Gyr.
    y=`python -c "print(0.2471 + ((0.2758-0.2471)*0.01/0.01817))"`
    echo "initial_y = $y" >> inlist_all
    echo "initial_mass = $w" >> inlist_all
    echo "0.01"
    echo '/ ! end of controls namelist' >> inlist_all
    srun --overlap -N 1 -o results/test-%j-$w-0.01.out -e results/test-%j-$w-0.01.err ./rn #Runs the script on the cluster with specified output and error file locations.
    wait
done

#Loop iterates over the masseshigh array
for w in "${masseshigh[@]}"
do 
    sleep 5
    echo "$w" #prints the current mass that is being simulated
    sed -i '$ d' inlist_high_mass #Removes the last 4 lines of the inlist_high_mass file (those corresponding to initial mass, helium abundance, and metallicity of the star)
    sed -i '$ d' inlist_high_mass
    sed -i '$ d' inlist_high_mass
    sed -i '$ d' inlist_high_mass
    echo "initial_z = 0.01" >> inlist_high_mass #Replaces those 4 lines with a star of mass w, metal content 0.01, and helium abundance that is solar scaled to 0.01 (Calculated by linearly scaling to the Helium abundance and metallicity that best reproduces solar observables at 4.61 Gyr.
    y=`python -c "print(0.2471 + ((0.2758-0.2471)*0.01/0.01817))"`
    echo "initial_y = $y" >> inlist_high_mass
    echo "initial_mass = $w" >> inlist_high_mass
    echo "0.01"
    echo '/ ! end of controls namelist' >> inlist_high_mass
    srun --overlap -N 1 -o results/test-%j-$w-0.01.out -e results/test-%j-$w-0.01.err ./rnhigh #Runs the ./rnhigh script on the cluster with specified output and error file locations.
    wait
done