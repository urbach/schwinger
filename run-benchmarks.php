<?php
	function callExecutableWithArguments($executable, $arguments)
	{
		$commandLine = "$executable $arguments > 'logs/$executable $arguments.log'";
		echo "Executing \"$commandLine\"...";
		$start = microtime(true);
		system($commandLine);
		$runtime = round(microtime(true) - $start, 3);
		echo " done (runtime: $runtime seconds).\n";
	}
	
	function buildArguments(&$argumentsList, $parameters, $parameterNames, $arguments)
	{
		if (!count($parameterNames))
		{
			array_push($argumentsList, trim($arguments));
			return;
		}
		
		$curParameterName = array_shift($parameterNames);
		
		foreach ($parameters[$curParameterName] as $curParameter)
			buildArguments($argumentsList, $parameters, $parameterNames, $arguments . ' --' . $curParameterName . ' ' . $curParameter);
	}
	
	$parameters = array(
			'beta' => array(
					'1.0 --mass -0.23125',
					/*'1.0 --mass -0.22750',
					'1.0 --mass -0.00600',
					'1.0 --mass -0.05750',
					'1.0 --mass -0.05500',
					'1.0 --mass 0.53750',
					'1.0 --mass 0.55000',
					'2.0 --mass -0.13250',
					'2.0 --mass -0.13125',
					'2.0 --mass -0.01250',
					'2.0 --mass -0.00750',
					'2.0 --mass 0.40000',
					'2.0 --mass 0.40625',
					'3.0 --mass -0.08250',
					'3.0 --mass -0.08188',
					'3.0 --mass -0.08125',
					'3.0 --mass 0.01875',
					'3.0 --mass 0.02000',
					'3.0 --mass 0.02125',
					'3.0 --mass 0.35000',
					'3.0 --mass 0.35625',
					'4.0 --mass -0.06000',
					'4.0 --mass -0.05600',
					'4.0 --mass 0.02750',
					'4.0 --mass 0.03000',
					'4.0 --mass 0.03125',
					'4.0 --mass 0.31250',
					'4.0 --mass 0.31500',
					'5.0 --mass -0.04375',
					'5.0 --mass -0.04250',
					'5.0 --mass 0.03250',
					'5.0 --mass 0.03500',
					'5.0 --mass 0.03750',
					'5.0 --mass 0.28500',
					'5.0 --mass 0.29000',
					'6.0 --mass -0.03505',
					'6.0 --mass -0.03250',
					'6.0 --mass 0.03000',
					'6.0 --mass 0.03750',
					'6.0 --mass 0.03775',
					'6.0 --mass 0.26000',
					'6.0 --mass 0.26250',
					'6.0 --mass 0.26500',*/
				),
			//'musqr' => array(0.1, 0.21, 0.4, 0.8),//array(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, '1.0'),
			//'n_steps_2' => array(2, 3, 4, 8, 16),//array(2, 4, 8, 12, 16),
			//'n_steps_1' => array(1, 2, 4),//array(1, 2, 4, 8),
			//'n_steps_0' => array(8),//array(1, 2, 4, 8),
		);
	
	$argumentsList = array();
	buildArguments($argumentsList, $parameters, array_keys($parameters), '');
	
	echo 'Total number of commands to execute: ' . count($argumentsList) . "\n";
	
	putenv('PATH=' . getenv('PATH'). ':~/Library/Developer/Xcode/DerivedData/qed-blwhrkwjolhqasfcwocujlitrtcl/Build/Products/Release');
	
	array_shift($argv);
	$defaultCommandLine = implode($argv, ' ');
	foreach ($argumentsList as $curArguments)
		callExecutableWithArguments('qed', $curArguments . ' ' . $defaultCommandLine);
?>