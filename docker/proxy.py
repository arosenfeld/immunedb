import argparse
import mimetypes
import os

import bottle

from immunedb.api.rest_service import create_app, EnableCors

mounted = set()
app = bottle.default_app()
app.install(EnableCors())


def mount_api(db_name):
    if db_name in mounted:
        return

    perms = {
        'database': db_name,
        'host': 'localhost',
        'username': 'root',
        'password': ''
    }

    app.mount(
        '/api/' + db_name + '/',
        create_app(perms)
    )
    mounted.add(db_name)


def mount_app(*args, **kwargs):
    bottle.request.environ['PATH_INFO'] = bottle.request.environ[
        'PATH_INFO'].rstrip('/')
    subdir = bottle.request.environ.get('HTTP_X_SCRIPT_NAME', '')
    if bottle.request.environ['PATH_INFO'].startswith(subdir):
        bottle.request.environ['PATH_INFO'] = (
            bottle.request.environ['PATH_INFO'][len(subdir) + 1:])

    path = bottle.request.environ['PATH_INFO'].strip('/').split('/')
    bottle.request.environ['PATH_INFO'] = os.path.join('/', *path)
    if len(path) < 2:
        bottle.abort(404)

    prefix, db_name = path[:2]
    if prefix == 'api':
        mount_api(db_name)


def replace_tokens(text, db_name):
    base = bottle.request.environ.get('HTTP_X_SCRIPT_NAME', '')[1:]
    replacements = {
        '__BASENAME__': os.path.join(base, 'frontend', db_name),
        '__ENDPOINT__': os.path.join('/' + base, 'api', db_name)
    }
    for search, replace in replacements.items():
        text = text.replace(search, replace)
    return text


@app.route('/frontend/<db_name>')
@app.route('/frontend/<db_name>/<filepath:path>')
def frontend(db_name, filepath=None):
    if not filepath:
        return bottle.redirect(
            os.path.join(
                '/frontend',
                db_name,
                'samples'))

    root = app.config['FRONTEND_ROOT']
    fn = os.path.join(root, filepath)
    if not os.path.exists(fn):
        fn = os.path.join(root, 'index.html')
    mime = mimetypes.guess_type(fn)[0]
    if mime and (mime == 'application/javascript' or mime.startswith('text')):
        with open(fn, encoding='utf-8') as fh:
            return replace_tokens(fh.read(), db_name)
    return bottle.static_file(filepath, root=root)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--frontend-path', help='Path to frontend',
                        default='/apps/immunedb-frontend')
    args = parser.parse_args()

    path = os.path.join(args.frontend_path, 'dist')
    if not os.path.exists(path):
        print('Path does not exist: {}'.format(path))
    else:
        app.add_hook('before_request', mount_app)
        app.config['FRONTEND_ROOT'] = path
        app.run(host='0.0.0.0', port=8080, debug=True)
