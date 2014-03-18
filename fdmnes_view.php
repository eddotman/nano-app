<?php
$handle = fopen("xanout/test_stand/cu_conv.txt", "r");

$res = array();


if ($handle) {
    while (($line = fgets($handle)) !== false) {
        // process the line read.
    	$nums = preg_split('/\s+/', $line);
    	array_shift($nums);
    	array_pop($nums);
    	array_push($res, array("x" => $nums[0] , "y" => $nums[1]));
    }
} else {
    // error opening the file.
}

array_shift($res);

echo json_encode($res);
?>