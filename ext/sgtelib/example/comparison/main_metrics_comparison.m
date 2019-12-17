close all
clear all
clc


output_dir = '.';



input_file = 'output_hartman3.txt';
output_file = [output_dir 'metric_comparison_hartman3.pdf'];
plot_metrics_comparison(input_file,output_file);
close all;

input_file = 'output_hartman6.txt';
output_file = [output_dir 'metric_comparison_hartman6.pdf'];
plot_metrics_comparison(input_file,output_file);
close all;

input_file = 'output_braninhoo.txt';
output_file = [output_dir 'metric_comparison_braninhoo.pdf'];
plot_metrics_comparison(input_file,output_file);
close all;

input_file = 'output_camelback.txt';
output_file = [output_dir 'metric_comparison_camelback.pdf'];
plot_metrics_comparison(input_file,output_file);
close all;

input_file = 'output_rosenbrock.txt';
output_file = [output_dir 'metric_comparison_rosenbrock.pdf'];
plot_metrics_comparison(input_file,output_file);
close all;
