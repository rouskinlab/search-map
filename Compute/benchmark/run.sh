#!/bin/zsh

### Benchmark the accurary of SEISMIC-RNA ###

set -euo pipefail


declare -A nprofiles=( [1]=1 [2]=4 [3]=3 [4]=3 )

# Simulate structures for each reference sequence.
for length in 280 560 1120
do
	ref=ref-$length
	pdir=$PWD/sim/params/$ref/full
	declare -A param3=( ["frag1"]=$((50. / $length + 0.5)) ["frag2"]=$((100. / $length + 0.5)) ["ampl2"]=1.0 )
	declare -A paraml=( ["frag1"]=$((100. / $length)) ["frag2"]=$((200. / $length)) ["ampl2"]=1.0 )
	seismic +sim ref --reflen $length --refs $ref --ref $ref
	fastqs=()
	for clusters in 1 2 3 4
	do
		cname=c$clusters
		ct=$pdir/$cname.ct
		seismic +sim fold --profile-name $cname --fold-max $clusters sim/refs/$ref.fa
		for library in frag1 frag2 ampl2
		do
			lct=$pdir/$library.ct
			if [ ! -e $lct ]
			then
				ln -s $ct $lct
			fi
			ends=$pdir/$library.ends.csv
			if [ ! -f $ends ]
			then
				p3=${param3[$library]}
				pl=${paraml[$library]}
				seismic +sim ends -3 $p3 -l $pl -i $lct
			fi
			for m in 4
			do
				mname=$cname-m$m
				mct=$pdir/$mname.ct
				if [ ! -e $mct ]
				then
					ln -s $ct $mct
				fi
				muts=$pdir/$mname.muts.csv
				if [ ! -f $muts ]
				then
					unm="0.0$m"
					seismic +sim muts -u am $unm -u cm $unm -i $mct
				fi
				for ((profile=1; profile<=${nprofiles[$clusters]}; profile+=1))
				do
					pname="c$clusters-$profile-m$m-$library"
					pct=$pdir/$pname.ct
					if [ ! -e $pct ]
					then
						ln -s $ct $pct
					fi
					if [ ! -e $pdir/$pname.ends.csv ]
					then
						ln -s $ends $pdir/$pname.ends.csv
					fi
					if [ ! -e $pdir/$pname.muts.csv ]
					then
						ln -s $muts $pdir/$pname.muts.csv
					fi
					if [ ! -e $pdir/$pname.clusts.csv ]
					then
						ln -s $PWD/clusts/c$clusters-$profile.csv $pdir/$pname.clusts.csv
					fi
					for reads in 10000 30000 100000 300000
					do
						sample=$pname-n$reads
						echo "Sample $sample for ref $ref"
						if [ $library = frag1 ]
						then
							fq=sim/samples/$sample/${ref}.fq.gz
							fastqs+=(-Z $fq)
							if [ ! -f $fq ]
							then
								seismic -v +sim fastq -s $sample -d $pdir -P $pname -n $reads --read-length 151 --single-end
							fi
						else
							fq1=sim/samples/$sample/${ref}_R1.fq.gz
							fq2=sim/samples/$sample/${ref}_R2.fq.gz
							fastqs+=(-X $fq1 -X $fq2)
							if [ $library = ampl2 ]
							then
								readlength=$(($length / 2 + 11))
							else
								readlength=151
							fi
							if [ ! -f $fq1 ] || [ ! -f $fq2 ]
							then
								seismic -v +sim fastq -s $sample -d $pdir -P $pname -n $reads --read-length $readlength --paired-end
							fi
						fi
					done
				done
			done
		done
	done
	# Process the simulated FASTQ files with SEISMIC-RNA.
	echo "Running SEISMIC-RNA on ${fastqs[@]}"
	seismic wf ${fastqs[@]} --no-fastqc --bt2-X $length -k 6 --no-table-read --no-graph-mprof --no-graph-tmprof --no-graph-ncov --no-graph-mhist --no-graph-giniroll --no-html sim/refs/$ref.fa
done

