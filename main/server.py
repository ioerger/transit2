# encoding: utf-8
from aiohttp import web
import socketio 
import sys
import argparse 
import math
import asyncio
from time import time as now


# 
# args setup
# 
parser = argparse.ArgumentParser(description="aiohttp server") 
parser.add_argument('--port') 
args = parser.parse_args()

# globals 
debugging = False

# 
# server setup
# 
app = web.Application(client_max_size=(1024 ** 2 * 100))
routes = web.RouteTableDef()
options = {} if not debugging else dict(logger=True, engineio_logger=True)
sio = socketio.AsyncServer(async_mode='aiohttp', cors_allowed_origins="*", **options); sio.attach(app)
cards = "{}"
large_data = {}

if debugging: print('starting server') 

import json
from os.path import join, dirname
__dir__ = dirname(__file__)

@routes.get('/')
async def index(request : web.Request): 
    try:
        html_path       = join(__dir__, 'index.html')
        javascript_path = f"{__dir__}/bundle.js"
        with open(html_path, 'r') as f:
            html = f.read()
        with open(javascript_path, 'r') as f:
            javascript = f.read()
        output = html.replace("##REPLACEME##", javascript)
    except Exception as error:
        print(f'''error = {error}''')
        output = None
    return web.Response(text=output, content_type='text/html')

@routes.post('/read')
async def read(request):
    arg = None
    try:
        arg = json.loads(await request.text())
    except Exception as error:
        print('error = ', error)
    
    try:
        with open(arg,'r') as f:
            output = f.read()
    except:
        output = None
    return web.Response(text=output)


@routes.post('/large/set/{content_type}/{data_id}')
async def set_large_data(request : web.Request):
    global large_data
    content_type = request.match_info["content_type"]
    content_type = content_type.replace(r"%2F", "/")
    large_data_id = request.match_info["data_id"]
    # save in ram
    post_result = await request.post()
    large_file = post_result.get("file")
    if large_file is not None:
        large_data[large_data_id] = large_file.file.read()
    return web.Response(text="null")

@routes.get('/large/get/{content_type}/{data_id}')
async def get_large_data(request : web.Request):
    global large_data
    content_type = request.match_info["content_type"]
    content_type = content_type.replace(r"%2F", "/")
    large_data_id = request.match_info["data_id"]
    return web.Response(
        content_type=content_type,
        body=large_data[large_data_id],
    )

# 
# start server
# 
app.add_routes(routes); web.run_app(app, port=args.port) 