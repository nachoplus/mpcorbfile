from configparser import ConfigParser
import sh

from pathlib import Path

# if pypirc is at the same directory level
PYPIRC_FILE = f"{Path.home()}/.pypirc"
print(f"getting credentials from: {PYPIRC_FILE}")


def main():
    parser = ConfigParser()
    parser.read(PYPIRC_FILE)

    # parse the section info
    pypirc = parser._sections
    index_servers = (pypirc.get("distutils")).get("index-servers").split()

    for server in index_servers:
        server_config = pypirc.get(server)
        print(f"Server:{server}")
        if server_config.get("repository"):
            sh.poetry.config(f"repositories.{server}", server_config["repository"])

        if server_config.get("username") and server_config.get("password"):
            print(f"Username:{server_config.get('username')}")
            sh.poetry.config(
                f"http-basic.{server}",
                server_config.get("username"),
                server_config.get("password"),
            )


if __name__ == "__main__":
    main()
