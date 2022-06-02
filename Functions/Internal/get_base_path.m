function base_path = get_base_path()
%GET_BASE_PATH returns the top level path environment variable
base_path = getenv("BASE_PATH");

if isempty(base_path)
    error("BASE_PATH not set, ensure setup function has been run.")
end
end