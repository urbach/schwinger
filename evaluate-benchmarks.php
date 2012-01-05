<?php
	function fitData($filename, $pion_measurements, $pion_errors, $pcac_measurements, $pcac_errors, $beta, $mass, $n_steps, $musqr)
	{
		$descriptorspec = array(
							   0 => array("pipe", "r"),  // STDIN ist eine Pipe, von der das Child liest
							   1 => array("pipe", "w"),  // STDOUT ist eine Pipe, in die das Child schreibt
							   2 => array("pipe", "w"),  // STDERR ist eine Pipe, in die das Child schreibt
		);
		
		$process = proc_open('gnuplot', $descriptorspec, $pipes, NULL, NULL);
		$gnuplot_input = $pipes[0];
		$gnuplot_output = $pipes[1];
		$gnuplot_err = $pipes[2];
		
		$nMeasurements = count($pion_measurements);
		
		$fitc = 0.5 * ($nMeasurements);// - 1);
		$fitstart = 0.25 * $nMeasurements;
		$fitend = 0.75 * $nMeasurements;
		
		for ($i = 0; $i < $nMeasurements; $i ++)
		{
			$pcac_diff[$i] = 0.5 * ($pcac_measurements[($i+1)%$nMeasurements] + $pcac_measurements[($i-1+$nMeasurements)%$nMeasurements]);
			
			$data_pion .= "$i\t$pion_measurements[$i]\t$pion_errors[$i]\n";
			$data_pcac .= "$i\t$pcac_measurements[$i]\t$pcac_errors[$i]\n";
			$data_pion_pcac .= "$i\t" . ($pion_measurements[$i] != 0 ? $pcac_measurements[$i] / (2 * $pion_measurements[$i]) : NAN) . "\n";
		}
		$fermion_mass = 0;
		for ($i = $fitstart; $i <= $fitend; $i ++)
			$fermion_mass += ($pion_measurements[$i] != 0 ? $pcac_measurements[$i] / (2 * $pion_measurements[$i]) : NAN);
		//echo $fermion_mass / ($fitend - $fitstart + 1) . "\n";
		$data_pion .= 'e';
		$data_pcac .= 'e';
		$data_pion_pcac .= 'e';
		for ($i = 0; $i < $nMeasurements - 1; $i ++)
		{
			$m1 = $pion_measurements[$i + 1];
			$m2 = $pion_measurements[$i];
			$data_pion_relation .= "$i\t" . ($m1 != 0 ? $m2 / $m1 : NAN) . "\t" . (($m1 != 0 && $m2 != 0) ? sqrt(($pion_errors[$i] / $m1)^2 + ($pion_errors[$i + 1] * $m2 / ($m1^2))^2) : NAN) . "\n";
		}
		$data_pion_relation .= 'e';
		$commands = <<<END
set font "HelveticaNeue"
set title "Mean correlation function over time, {/Symbol b} = $beta, m = $mass, n_{steps} = [$n_steps[0], $n_steps[1], $n_steps[2]], {/Symbol m}^2 = $musqr"
set xlabel "t"
set ylabel "C(t)"
set xrange [0:$nMeasurements]
set terminal postscript enhanced "HelveticaNeue" color eps
set output "plots/$filename.eps"
#set bars small
set logscale y
b=0.2
c=$fitc
f(x)=0.5*a*(exp(-b*x)+exp(-b*(2*c-x)))
fit [$fitstart:$fitend] f(x) '-' using 1:2 via a, b
$data_pion
plot '-' using 1:2:3 with yerrorbars title "C_{/Symbol p}(t)", f(x) title "Fit" with lines, '-' using 1:2:3 with yerrorbars title "C(t) / C(t + 1)", '-' using 1:2:3 with yerrorbars title "C_{PCAC}(t)", '-' using 1:2 title "C_{PCAC}(t) / (2 * C_{/Symbol p}(t))"
$data_pion
$data_pion_relation
$data_pcac
$data_pion_pcac
END;
		
		fprintf($gnuplot_input, $commands);
		fclose($gnuplot_input);
		while (stristr($gnuplot_results, 'correlation matrix of the fit parameters:') === FALSE)
			$gnuplot_results .= fread($gnuplot_err, 1000000);
		fclose($gnuplot_err);
		fclose($gnuplot_output);
		
		$results = explode('Final set of parameters', $gnuplot_results);
		preg_match('/b \s*= ([0-9\.]+)\s*\+\/- ([0-9\.]+)\s*\(/', $results[1], $matches);
		array_shift($matches);
		return $matches;
	}

	function extractDataFromLog($logname)
	{
		$logtext = file_get_contents($logname);
		$parts = explode(' Algorithm configuration:', $logtext);
		$results = $parts[1];
		
		// parameters
		preg_match('/no_timescales = (.)/', $parts[0], $matches);
		if (count($matches))
			$no_timescales = $matches[1];
		else
			$no_timescales = 3;
		
		preg_match('/X_1 = (.+), beta = (.+), g_mass = (.+), g_musqr = (.+)/', $results, $matches);
		array_shift($matches);
		list($X1, $beta, $mass, $musqr) = $matches;
		
		preg_match('/g_measurements = (.+), g_intermediate = (.+), n_steps\[0\] = (.+), n_steps\[1\] = (.+), n_steps\[2\] = (.+)/', $results, $matches);
		array_shift($matches);
		list($measurements, $intermediate, $n_steps[0], $n_steps[1], $n_steps[2]) = $matches;
		
		
		// performance
		preg_match('/Acceptance rate:\s*(.+)/', $results, $matches);
		array_shift($matches);
		list($acceptance_rate) = $matches;
		
		preg_match_all('/CG iterations per.*:\s*(.+)/', $results, $matches);
		$matches = $matches[1];
		list($cg_inner_update, $cg_outer_update, $cg_inner_solve, $cg_outer_solve) = $matches;
		
		preg_match_all('/force per.*:\s*([0-9\.]+)/', $results, $matches);
		$matches = $matches[1];
		list($force_gauge_update, $force_inner_update, $force_outer_update, $force_gauge_application, $force_inner_application, $force_outer_application) = $matches;
		
		preg_match('/Runtime \/ seconds:\s*(.+)/', $results, $matches);
		array_shift($matches);
		list($runtime) = $matches;
		
		
		// measurements
		preg_match('/Plaquette:\s*(.+) \+\/- (.+) \(/', $results, $matches);
		array_shift($matches);
		list($plaquette, $plaquette_e) = $matches;
		
		preg_match('/Plaquette autocorrelation time:\s*(.+)/', $results, $matches);
		array_shift($matches);
		list($plaquette_autocorrelation_time) = $matches;
		
		preg_match('/Polyakov loop:\s*(.+) \+\/- (.+) \(/', $results, $matches);
		array_shift($matches);
		list($pl, $pl_e) = $matches;
		
		preg_match('/Chiral Condensate:\s*(.+) \+\/- (.+) \(/', $results, $matches);
		array_shift($matches);
		list($cc, $cc_e) = $matches;
		
		preg_match('/-Delta H:\s*(.+) \+\/- (.+) \(/', $results, $matches);
		array_shift($matches);
		list($dh, $dh_e) = $matches;
		
		preg_match('/exp(-Delta H):\s*(.+) \+\/- (.+) \(/', $results, $matches);
		array_shift($matches);
		list($expdh, $expdh_e) = $matches;
		
		preg_match_all('/Pion Correlation\[(.+)\]:\s*(.+) \+\/- (.+) \(/', $results, $matches);
		array_shift($matches);
		array_shift($matches);
		list($pion_correlation, $pion_correlation_e) = $matches;
		
		preg_match_all('/PCAC Correlation\[(.+)\]:\s*(.+) \+\/- (.+) \(/', $results, $matches);
		array_shift($matches);
		array_shift($matches);
		list($pcac_correlation, $pcac_correlation_e) = $matches;
		
		if (count($pion_correlation))
			list($m_pion, $m_pion_e) = fitData(basename($logname), $pion_correlation, $pion_correlation_e, $pcac_correlation, $pcac_correlation_e, $beta, $mass, $n_steps, $musqr);
		
		return array(basename($logname), $X1, $beta, $mass,
					 $no_timescales, $n_steps[0], $n_steps[1], $n_steps[2], $musqr,
					 
					 $acceptance_rate, $runtime,
					 $cg_inner_update, $cg_outer_update,
					 $force_inner_application, $force_outer_application,
					 
					 $plaquette, $plaquette_autocorrelation_time, $pl, $cc,
					 $m_pion, $m_pion_e);
	}
	
	array_shift($argv);
	
	if (count($argv))
		$lognames = $argv;
	else
		$lognames = explode("\n", file_get_contents('filelist.txt'));
	
	echo implode("\t", array('filename', 'X_1', 'beta    ', 'mass    ',
							 'scales', 'steps_0', 'steps_1', 'steps_2', 'musqr    ',

							 'acc', 'runtime',
							 'cg_in', 'cg_out',
							 'force_inner', 'force_outer',

							 'plaquette', 'plaquette_time', 'polyakov', 'cc',
							 'm_pion', 'm_pion_e')) . "\n";
	
	foreach ($lognames as $logname)
		if (strlen($logname) && substr($logname, 0, 1) != "#")
			echo implode("\t", extractDataFromLog($logname)) . "\n";
?>