import subprocess




# interface
class Executor:


    # function to run a cmd
    def runCmd (self, cmdString:str) -> None:
        assert (False)


# implementation of Executor class that runs the commands directly
class SubprocessExecutor (Executor):


    def __init__(self) -> None:
        # nothing to do
        pass


    def __del__ (self) -> None:
        # nothing to do
        pass


    # inherited: function to run a cmd
    def runCmd (self, cmdString:str) -> None:
        # use subprocess to run the command
        print (cmdString)
        result = subprocess.run (cmdString, shell=True, capture_output=True, text=True)
        assert ((result.returncode == 0) and (result.stderr == '' or result.stderr is None)), (result.returncode, result.stderr)
        print (result.stdout)


# implementation of Executor class that appends all the commands to a file
class FileExecutor (Executor):


    def __init__(self, filename:str, append:bool=False) -> None:
        # to append or not to append
        if (append):
            self._fileMode = 'a'
        else:
            self._fileMode = 'w'
        self._filename = filename
        self._ofs = None

    def __del__ (self) -> None:
        # close the file if it was ever opened
        if (self._ofs is not None):
            self._ofs.close()


    # inherited: function to run a cmd
    def runCmd (self, cmdString:str) -> None:

        # lazy evaluation: only open file if we write commands
        if (self._ofs is None):
            self._ofs = open (self._filename, self._fileMode)
            # better save then sorry
            if (self._fileMode == 'a'):
                self._ofs.write ('\n')

        # write the command into file
        self._ofs.write (cmdString + '\n')

