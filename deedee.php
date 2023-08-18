<?php
session_start();

$email = $_POST["email"];
$assay = $_POST["assay"];
$gene = $_POST["gene"];

if(isset($_SESSION['upD'])){
   $upload_dir = $_SESSION['upD'];
   $run_comand = "sbatch --ntasks=1 --cpus-per-task 2 --error " . escapeshellarg($upload_dir) . "/run_log.e --output " . escapeshellarg($upload_dir) . "/run_log.o /var/www/html/scripts/deedee/deedee.sh " . escapeshellarg($upload_dir) . " " . $email . " " . $assay . " " . $gene;
   shell_exec($run_comand);
   unset($_SESSION['upD']);
   session_regenerate_id(true);
   header("Location: /uploaded.html");
  } else{
     session_regenerate_id(true);
     header("Location: /NOTuploaded.html");
  }
?>