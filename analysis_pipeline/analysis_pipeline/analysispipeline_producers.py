import os, sys
from analytics_utils.database_access.db_interface import DatabaseInterface
from analytics_utils.utils.slack_interface import post_slack_message, SlackMessageConstants
import json


class Environment:
    DEV = 'DEV'
    PROD = 'PROD'

    
class MessageProducer:
    ML = 'ML'
    LAB_AUTO = 'LAB_AUTO'
    

class BifröstMessageProducer:
    def __init__(self, sender, env):
        self.sender = sender
        self.env = env
        
    def add_event(self, payload, queue_name):
        post_slack_message(json.dumps(payload), queue_name)        


ML_producer = BifröstMessageProducer(sender=MessageProducer.ML, env=Environment.PROD)
LabAuto_producer = BifröstMessageProducer(sender=MessageProducer.LAB_AUTO, env=Environment.PROD)