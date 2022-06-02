function results_path = get_results_path()
%GET_RESULTS_PATH returns the top level path environment variable
results_path = getenv("RESULTS_PATH");

if isempty(results_path)
    error("RESULTS_PATH not set, ensure setup function has been run.")
end
end

