from jumpssh import SSHSession

gateway_session = SSHSession('styx.jb.man.ac.uk', 'gsearle', password='').open()
remote_session = gateway_session.get_remote_session('slioch', password='')
print(remote_session.get_cmd_output('ls'))
#remote_session.get('/remote/path/to/the/file', '/local/path/')
print(range(1,9))
