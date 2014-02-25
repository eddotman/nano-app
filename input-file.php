<?php
	$title = $_POST['title'];
	$outfile = $_POST['outfile'];
	$range = $_POST['range'];
	$radius = $_POST['radius'];
	$crystal = $_POST['crystal'];
	$atom = $_POST['atom'];
	$efermi = $_POST['efermi'];

	$arr = array(
		'title_line' => $title,
		'outfile' => $outfile,
		'range' => $range,
		'radius' => $radius,
		'crystal_dim' => $crystal,
		'atoms' => $atom,
		'fermi_energy' => $efermi
		);

	$json = json_encode($arr);

	$output = array();
	system("python create-input-file.py '" . $json . "'");

	echo "inp/" . $title . ".inp";
?>